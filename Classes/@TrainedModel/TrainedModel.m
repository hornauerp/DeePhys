classdef TrainedModel
% TRAINEDMODEL  Serializable deployment wrapper for a trained classifier/regressor.
%
%   Bundles a trained model with its normalization pipeline and feature
%   specification, so new data can be processed through the exact same
%   pipeline in one call.
%
%   USAGE:
%     % After training:
%     tm = TrainedModel.fromClassificationResult(result, np, rg);
%     save('my_model.mat', 'tm');
%
%     % Later, on new data:
%     tm = load('my_model.mat').tm;
%     [Y_pred, scores] = tm.predict(new_recordings);
%
%   The object stores everything needed to reproduce the prediction pipeline:
%   feature groups, normalization parameters, and the classifier itself.

    properties
        Mdl                     % Trained MATLAB classifier/regressor (compact or full)
        Normalization           % NormalizationPipeline with fitted FitParams
        FeatureSpec             % Struct: .level, .unit_features, .network_features, .feature_names, .useClustered
        TaskType     string     % "classification" or "regression"
        GroupLabels             % Class labels (classification) or regression variable name
        TrainDate    datetime   % When the model was trained
        Parameters   struct     % All training parameters for provenance
        Preprocessing           % Fold preprocessing from MLResult (medians, pipeline, column mask);
                                % when set, predict() replays it instead of Normalization
    end

    methods

        function tm = TrainedModel(mdl, normalization, feature_spec, task_type, group_labels, parameters)
        % Constructor.
            arguments
                mdl = []
                normalization = []
                feature_spec struct = struct()
                task_type string = "classification"
                group_labels = []
                parameters struct = struct()
            end
            tm.Mdl = mdl;
            tm.Normalization = normalization;
            tm.FeatureSpec = feature_spec;
            tm.TaskType = task_type;
            tm.GroupLabels = group_labels;
            tm.TrainDate = datetime('now');
            tm.Parameters = parameters;
        end

        function [Y_pred, scores] = predict(tm, new_objects, normalization_var, group_labels_for_norm)
        % PREDICT  Apply the full pipeline to new objects and return predictions.
        %
        % INPUTS:
        %   new_objects             - FeatureStore, RecordingProcessor array, feature table
        %                             (columns named as in training), or the legacy
        %                             MEArecording / Unit arrays
        %   normalization_var       - (optional) metadata field for group labels in normalization.
        %                             Only needed if the pipeline uses group normalization.
        %   group_labels_for_norm   - (optional) explicit group labels for new_objects.
        %                             If not provided, derived from metadata.
        %
        % OUTPUTS:
        %   Y_pred - predicted labels (classification) or values (regression)
        %   scores - prediction scores/probabilities (classification) or [] (regression)
            arguments
                tm TrainedModel
                new_objects
                normalization_var string = ""
                group_labels_for_norm = []
            end

            % Step 1: Extract features using the same specification
            fs = tm.FeatureSpec;
            if istable(new_objects)
                input_table = new_objects;
            elseif isa(new_objects, 'FeatureStore') || isa(new_objects, 'RecordingProcessor')
                fstore = new_objects;
                if isa(fstore, 'RecordingProcessor')
                    fstore = FeatureStore.fromProcessors(fstore);
                end
                if isfield(fs, 'level') && fs.level == "Unit"
                    input_table = fstore.unitMatrix(TrainedModel.specGroups(fs, 'unit_features'));
                else
                    input_table = fstore.recordingMatrix(TrainedModel.specGroups(fs, 'network_features'));
                end
            elseif isa(new_objects, 'Unit')
                input_table = FeatureAssembly.unitFeatures(new_objects, fs.unit_features);
            elseif isa(new_objects, 'MEArecording')
                if fs.level == "Unit"
                    input_table = FeatureAssembly.unitFeatures(new_objects, fs.unit_features);
                else
                    input_table = FeatureAssembly.recordingFeatures(new_objects, fs.network_features, fs.unit_features, fs.useClustered);
                end
            else
                error('TrainedModel:predict', ['Input must be a FeatureStore, RecordingProcessor ' ...
                    'array, feature table, or a legacy MEArecording / Unit array.']);
            end

            if ~isempty(tm.Preprocessing)
                % Step 2: replay the training fold's imputation, normalisation and column mask
                % (columns are matched by name, so input column order does not matter)
                X = MLUtils.applyPreprocessing(tm.Preprocessing, input_table, group_labels_for_norm);
            else
                % Select features by name, in training order
                if ~isempty(fs.feature_names)
                    vars = string(input_table.Properties.VariableNames);
                    [found, loc] = ismember(string(fs.feature_names), vars);
                    if ~all(found)
                        error('TrainedModel:predict', 'New data lacks %d training feature(s).', sum(~found));
                    end
                    input_table = input_table(:, loc);
                end

                % Step 2: Normalize using the fitted pipeline
                X = input_table.Variables;
                X(isnan(X)) = 0;

                if ~isempty(tm.Normalization) && ~isempty(fieldnames(tm.Normalization.FitParams))
                    X = tm.Normalization.transform(X, group_labels_for_norm);
                end
            end

            % Step 3: Predict
            if tm.TaskType == "classification"
                [Y_pred, scores] = tm.Mdl.predict(X);
                if ~isempty(tm.Preprocessing) && ~isempty(tm.GroupLabels) && isnumeric(Y_pred)
                    Y_pred = tm.GroupLabels(Y_pred);   % Classifier.classify trains on class indices
                end
            else
                Y_pred = tm.Mdl.predict(X);
                scores = [];
            end
        end

        function summary = describe(tm)
        % DESCRIBE  Return a human-readable summary of the trained model.
            summary = sprintf('TrainedModel (%s)\n  Trained: %s\n  Features: %s\n  Level: %s\n', ...
                tm.TaskType, string(tm.TrainDate), ...
                strjoin(string(fieldnames(tm.FeatureSpec)), ', '), ...
                tm.FeatureSpec.level);
        end

    end

    methods (Static)

        function tm = fromClassificationResult(clf_result, normalization_pipeline, parameters)
        % FROMCLASSIFICATIONRESULT  Build a TrainedModel from a ClassificationResult.
        %
        % INPUTS:
        %   clf_result              - ClassificationResult (single fold, typically best or final)
        %   normalization_pipeline  - fitted NormalizationPipeline
        %   parameters              - struct with training parameters (level, features, etc.)
            arguments
                clf_result ClassificationResult
                normalization_pipeline NormalizationPipeline = NormalizationPipeline()
                parameters struct = struct()
            end
            feature_spec = struct();
            feature_spec.level = "Recording";
            feature_spec.unit_features = "all";
            feature_spec.network_features = "all";
            feature_spec.useClustered = false;
            feature_spec.feature_names = clf_result.Features;
            % Override from parameters if provided
            if isfield(parameters, 'level'),            feature_spec.level = parameters.level; end
            if isfield(parameters, 'unit_features'),    feature_spec.unit_features = parameters.unit_features; end
            if isfield(parameters, 'network_features'), feature_spec.network_features = parameters.network_features; end
            if isfield(parameters, 'useClustered'),     feature_spec.useClustered = parameters.useClustered; end

            tm = TrainedModel(clf_result.Mdl, normalization_pipeline, feature_spec, ...
                "classification", clf_result.GroupLabels, clf_result.Parameters);
            tm = TrainedModel.attachPreprocessing(tm, clf_result);
        end

        function tm = fromRegressionResult(reg_result, normalization_pipeline, parameters)
        % FROMREGRESSIONRESULT  Build a TrainedModel from a RegressionResult.
            arguments
                reg_result RegressionResult
                normalization_pipeline NormalizationPipeline = NormalizationPipeline()
                parameters struct = struct()
            end
            feature_spec = struct();
            feature_spec.level = "Recording";
            feature_spec.unit_features = "all";
            feature_spec.network_features = "all";
            feature_spec.useClustered = false;
            feature_spec.feature_names = reg_result.Features;
            if isfield(parameters, 'level'),            feature_spec.level = parameters.level; end
            if isfield(parameters, 'unit_features'),    feature_spec.unit_features = parameters.unit_features; end
            if isfield(parameters, 'network_features'), feature_spec.network_features = parameters.network_features; end
            if isfield(parameters, 'useClustered'),     feature_spec.useClustered = parameters.useClustered; end

            tm = TrainedModel(reg_result.Mdl, normalization_pipeline, feature_spec, ...
                "regression", [], reg_result.Parameters);
            tm = TrainedModel.attachPreprocessing(tm, reg_result);
        end

    end

    methods (Static, Access = private)

        function g = specGroups(spec, field)
        % Feature groups stored in the FeatureSpec ("all" when not recorded).
            g = "all";
            if isfield(spec, field) && ~isempty(spec.(field))
                g = string(spec.(field));
            end
        end

        function tm = attachPreprocessing(tm, ml_result)
        % Use the fold's own fitted preprocessing when the result carries it. The feature
        % list then becomes the pre-mask input names, which applyPreprocessing expects.
            if ~isempty(ml_result.Preprocessing)
                tm.Preprocessing = ml_result.Preprocessing;
                tm.FeatureSpec.feature_names = ml_result.Preprocessing.InputFeatures;
            end
        end

    end
end
