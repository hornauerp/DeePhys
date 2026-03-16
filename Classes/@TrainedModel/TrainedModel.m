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
        %   new_objects             - MEArecording array or Unit array
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
            if isa(new_objects, 'Unit')
                input_table = FeatureAssembly.unitFeatures(new_objects, fs.unit_features);
            elseif isa(new_objects, 'MEArecording')
                if fs.level == "Unit"
                    input_table = FeatureAssembly.unitFeatures(new_objects, fs.unit_features);
                else
                    input_table = FeatureAssembly.recordingFeatures(new_objects, fs.network_features, fs.unit_features, fs.useClustered);
                end
            else
                error('TrainedModel:predict', 'Input must be Unit or MEArecording array.');
            end

            % Select features by name if specified
            if ~isempty(fs.feature_names)
                vars = string(input_table.Properties.VariableNames);
                sel_vars = ismember(vars, fs.feature_names);
                input_table = input_table(:, sel_vars);
            end

            % Step 2: Normalize using the fitted pipeline
            X = input_table.Variables;
            X(isnan(X)) = 0;

            if ~isempty(tm.Normalization) && ~isempty(fieldnames(tm.Normalization.FitParams))
                if isempty(group_labels_for_norm)
                    group_labels_for_norm = [];
                end
                X = tm.Normalization.transform(X, group_labels_for_norm);
            end

            % Step 3: Predict
            if tm.TaskType == "classification"
                [Y_pred, scores] = tm.Mdl.predict(X);
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
        end

    end
end
