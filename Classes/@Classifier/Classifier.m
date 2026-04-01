classdef Classifier
% CLASSIFIER  Static classification methods that operate on plain feature tables.
%
% Zero coupling to MEArecording, Unit, or RecordingGroup. Accepts a feature table
% (from FeatureStore.unitMatrix / recordingMatrix / cultureMatrix) and a label
% vector Y directly.
%
% USAGE:
%   opts.Algorithm         = 'rf';          % 'rf', 'svm', 'cnb', 'knn'
%   opts.KFold             = -1;            % -1 = LOOCV
%   opts.NHyper            = 0;             % 0 = no hyperparam search
%   opts.CVGroups          = fs.UnitTable.RecordingID;   % recording-level CV grouping
%   opts.NormalizationGroups = fs.UnitTable.PlatingDate; % per-plate normalization
%   opts.PoolingValues     = {};            % merge class labels
%
%   result = Classifier.classify(X, Y, opts);
%   result = Classifier.classifyByGroups(X, Y, feature_group_defs, opts);
%   [Y_pred, scores] = Classifier.predict(trained_model, X);

    methods (Static)

        function result = classify(X, Y, opts)
        % CLASSIFY  K-fold cross-validated classification on a feature table.
        %
        % INPUTS:
        %   X    - (N x F) table of feature values (from fs.unitMatrix etc.)
        %   Y    - (N x 1) categorical / string / double label vector
        %   opts - struct with optional fields:
        %     .Algorithm         - 'rf' (default), 'svm', 'cnb', 'knn'
        %     .KFold             - number of folds; -1 = LOOCV (default -1)
        %     .NHyper            - hyperparam search evaluations (default 0)
        %     .CVGroups          - (N x 1) string/categorical for hierarchy-aware CV
        %                          (e.g. fs.UnitTable.RecordingID). [] = standard CV.
        %     .CVLevel           - 'recording' or 'culture' when CVGroups provided
        %     .NormalizationGroups - (N x 1) string/categorical for per-group z-score
        %                           (e.g. fs.UnitTable.PlatingDate). [] = global z-score.
        %     .PoolingValues     - cell array to merge label values
        %     .GroupLabels       - (optional) string array of class label names
        %
        % OUTPUTS:
        %   result - (1 x K) ClassificationResult array, one per fold
            arguments
                X    table
                Y
                opts struct = struct()
            end

            opts = Classifier.defaultOpts(opts);
            [Y, group_labels] = Classifier.prepareLabels(Y, opts.PoolingValues);

            norm_groups = Classifier.resolveGroups(opts.NormalizationGroups, height(X));
            cv_groups   = Classifier.resolveGroups(opts.CVGroups, height(X));

            % Build CV partition
            K = opts.KFold;
            if ~isempty(cv_groups)
                gcv     = GroupedCV.byGroups(cv_groups, K, Y);
                K       = gcv.NumFolds;
                use_gcv = true;
            else
                if K == -1
                    K = numel(Y);
                end
                cv      = cvpartition(Y, 'KFold', K);
                use_gcv = false;
            end

            t_start = tic;
            result(K) = ClassificationResult();
            for k = 1:K
                if use_gcv
                    [train_idx, test_idx] = gcv.fold(k);
                else
                    train_idx = cv.training(k);
                    test_idx  = cv.test(k);
                end

                Y_train = Y(train_idx);
                Y_test  = Y(test_idx);

                [X_train, X_test, feat_names] = Classifier.normalizeFeatures( ...
                    X, train_idx, test_idx, norm_groups);

                [clf, train_acc] = MLPipeline.createClassifier(X_train, Y_train, opts.Algorithm, opts.NHyper);
                [Y_pred, scores] = predict(clf, X_test);

                if opts.Algorithm == "rf"
                    predImp = oobPermutedPredictorImportance(clf, 'Options', statset('UseParallel', true));
                else
                    predImp = [];
                end

                result(k) = ClassificationResult( ...
                    'Mdl',         clf, ...
                    'Y_pred',      Y_pred, ...
                    'Y_test',      Y_test, ...
                    'scores',      scores, ...
                    'objects',     [], ...
                    'train_acc',   train_acc, ...
                    'predImp',     predImp, ...
                    'GroupLabels', group_labels, ...
                    'Features',    feat_names, ...
                    'Parameters',  struct( ...
                        'algorithm',        opts.Algorithm, ...
                        'normalization_var', opts.NormalizationGroups, ...
                        'K_fold',           K, ...
                        'N_hyper',          opts.NHyper));
            end

            elapsed = toc(t_start);
            summary = ClassificationResult.summarizeFolds(result);
            AnalysisLog.instance().add('Classifier.classify', opts, ...
                sprintf('Acc=%.3f±%.3f F1=%.3f', summary.mean_accuracy, summary.std_accuracy, summary.mean_F1), ...
                elapsed);
        end

        function results = classifyByGroups(X, Y, feature_group_defs, opts)
        % CLASSIFYBYGROUPS  Run classify() over a set of feature-column subsets.
        %
        % feature_group_defs is a struct array where each element has:
        %   .Name    - string label for the group
        %   .Columns - string array of column names (or a FeatureStore feature group name)
        %
        % Returns a struct array with .Name and .Result fields.
            arguments
                X                 table
                Y
                feature_group_defs struct
                opts               struct = struct()
            end
            results = struct('Name', {}, 'Result', {});
            for g = 1:numel(feature_group_defs)
                fg   = feature_group_defs(g);
                cols = string(fg.Columns);
                all_cols = string(X.Properties.VariableNames);
                sel = ismember(all_cols, cols) | startsWith(all_cols, cols);
                if ~any(sel)
                    warning('Classifier:classifyByGroups', ...
                        'No columns matched for group "%s" — skipped.', fg.Name);
                    continue
                end
                Xsub = X(:, sel);
                res  = Classifier.classify(Xsub, Y, opts);
                results(end+1) = struct('Name', fg.Name, 'Result', res); %#ok<AGROW>
            end
        end

        function [Y_pred, scores] = predict(trained_model, X)
        % PREDICT  Apply a trained classifier (TrainedModel or raw clf) to new data.
        %
        %   trained_model - TrainedModel object (from classifyByFeatureGroups) or
        %                   a raw MATLAB classifier
        %   X             - (N x F) table of feature values
            if isa(trained_model, 'TrainedModel')
                [Y_pred, scores] = trained_model.predict(X);
            else
                [Y_pred, scores] = predict(trained_model, X.Variables);
            end
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Static, Access = private)

        function opts = defaultOpts(opts)
            if ~isfield(opts, 'Algorithm'),          opts.Algorithm = 'rf';       end
            if ~isfield(opts, 'KFold'),              opts.KFold = -1;             end
            if ~isfield(opts, 'NHyper'),             opts.NHyper = 0;             end
            if ~isfield(opts, 'CVGroups'),           opts.CVGroups = [];          end
            if ~isfield(opts, 'CVLevel'),            opts.CVLevel = 'recording';  end
            if ~isfield(opts, 'NormalizationGroups'), opts.NormalizationGroups = []; end
            if ~isfield(opts, 'PoolingValues'),      opts.PoolingValues = {};     end
        end

        function [Y_out, group_labels] = prepareLabels(Y, pooling_vals)
        % Convert Y to a consistent format and apply any label pooling.
            if ischar(Y) || (iscell(Y) && ischar(Y{1}))
                Y = string(Y);
            end
            if ~isempty(pooling_vals) && iscell(pooling_vals) && ~isempty(pooling_vals{1})
                [Y_out, group_labels] = MLPipeline.poolMetadataValues(Y, unique(Y,'stable'), pooling_vals);
            else
                [group_labels, ~, Y_out] = unique(Y, 'stable');
            end
        end

        function g = resolveGroups(groups_input, N)
        % Normalise the groups input to a (N x 1) double index vector, or [].
            if isempty(groups_input)
                g = [];
                return
            end
            if ischar(groups_input) || isstring(groups_input)
                % String array of group labels → numeric index
                [~, ~, g] = unique(string(groups_input(:)), 'stable');
            elseif isnumeric(groups_input) || iscategorical(groups_input)
                [~, ~, g] = unique(groups_input(:), 'stable');
            else
                g = [];
            end
            if numel(g) ~= N
                error('Classifier:resolveGroups', ...
                    'Group vector length (%d) must match number of observations (%d).', numel(g), N);
            end
        end

        function [X_train, X_test, feat_names] = normalizeFeatures(X, train_idx, test_idx, norm_groups)
        % Impute NaN, fit normalization on training set, apply to test set.
        % Mirrors prepareInputMatrix but with no RecordingGroup dependency.
            mat = X.Variables;
            feat_names = string(X.Properties.VariableNames);

            % Impute NaN per column using training-set median
            for col = 1:size(mat, 2)
                train_col = mat(train_idx, col);
                med = median(train_col(~isnan(train_col)));
                if isnan(med), med = 0; end
                mat(isnan(mat(:, col)), col) = med;
            end

            % Normalization pipeline
            if ~isempty(norm_groups)
                np = NormalizationPipeline.groupThenGlobal();
                g_train = norm_groups(train_idx);
            else
                np = NormalizationPipeline.globalOnly();
                g_train = [];
            end

            [X_train, np] = np.fit_transform(mat(train_idx, :), g_train);

            if any(test_idx)
                g_test = [];
                if ~isempty(norm_groups)
                    g_test = norm_groups(test_idx);
                end
                X_test = np.transform(mat(test_idx, :), g_test);
            else
                X_test = [];
            end

            % Drop zero-variance / all-NaN columns
            bad = any(isnan(X_train));
            if ~isempty(X_test)
                bad = bad | any(isnan(X_test));
                X_test(:, bad) = [];
            end
            X_train(:, bad) = [];
            feat_names(bad) = [];
        end

    end
end
