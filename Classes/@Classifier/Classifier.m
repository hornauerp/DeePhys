classdef Classifier
% CLASSIFIER  Static classification methods that operate on plain feature tables.
%
% Zero coupling to MEArecording, Unit, or RecordingGroup. Accepts a feature table
% (from FeatureStore.unitMatrix / recordingMatrix / cultureMatrix) and a label
% vector Y directly.
%
% USAGE:
%   opts.Algorithm         = 'rf';          % 'rf', 'svm', 'cnb', 'knn'
%   opts.KFold             = 5;             % -1 = leave-one-group-out
%   opts.NHyper            = 0;             % 0 = no hyperparam search
%   opts.CVGroups          = fs.UnitTable.RecordingID;   % recording-level CV grouping
%   opts.NormalizationGroups = fs.UnitTable.PlatingDate; % per-plate normalization
%   opts.PoolingValues     = {};            % merge class labels
%
%   result = Classifier.classify(X, Y, opts);
%   result = Classifier.classifyByGroups(X, Y, feature_group_defs, opts);
%   [Y_pred, scores] = Classifier.predict(result(k).Mdl, X);

    methods (Static)

        function result = classify(X, Y, opts)
        % CLASSIFY  K-fold cross-validated classification on a feature table.
        %
        % INPUTS:
        %   X    - (N x F) table of feature values (from fs.unitMatrix etc.)
        %   Y    - (N x 1) categorical / string / double label vector
        %   opts - struct with optional fields:
        %     .Algorithm         - 'rf' (default), 'svm', 'cnb', 'knn'
        %     .KFold             - number of folds; -1 = leave-one-group-out (default 5)
        %     .NHyper            - hyperparam search evaluations (default 0)
        %     .CVGroups          - (N x 1) string/categorical for hierarchy-aware CV
        %                          (e.g. fs.UnitTable.RecordingID). [] = standard CV.
        %     .CVLevel           - 'recording' or 'culture' when CVGroups provided
        %     .NormalizationGroups - (N x 1) string/categorical for per-group z-score
        %                           (e.g. fs.UnitTable.PlatingDate). [] = global z-score.
        %     .PoolingValues     - cell array to merge label values
        %     .GroupLabels       - (optional) string array of class label names
        %     .Prior             - 'empirical' (default) or 'uniform' (equal class weight)
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

            % Reproducibility: fix rng seed before building CV partition.
            if ~isempty(opts.Seed)
                rng(opts.Seed);
            end

            % Warn if class distribution is severely imbalanced (>3:1 ratio).
            Classifier.checkClassBalance(Y);

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

                clf_params = MLPipeline.returnDefaultParams();
                clf_params.RF.Prior = opts.Prior;
                [clf, train_acc] = MLPipeline.createClassifier(X_train, Y_train, opts.Algorithm, opts.NHyper, clf_params);
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
                        'normalization_var', opts.NormalizationVar, ...
                        'K_fold',           K, ...
                        'N_hyper',          opts.NHyper, ...
                        'seed',             opts.Seed));
            end

            elapsed = toc(t_start);
            summary = ClassificationResult.summarizeFolds(result);
            AnalysisLog.instance().add('Classifier.classify', opts, ...
                sprintf('Acc=%.3f±%.3f F1=%.3f', summary.mean_accuracy, summary.std_accuracy, summary.mean_F1), ...
                elapsed);
        end

        function [p_value, null_accs, obs_acc] = permutationTest(X, Y, opts, n_permutations)
        % PERMUTATIONTEST  Estimate classification significance via label permutation.
        %
        % Shuffles class labels n_permutations times and runs classify() each time,
        % building a null distribution of accuracies. Returns a one-sided p-value:
        %   P(null_accuracy >= observed_accuracy)
        %
        % INPUTS:
        %   X              - (N x F) feature table
        %   Y              - (N x 1) label vector
        %   opts           - same opts struct as Classifier.classify (reused for each run)
        %   n_permutations - number of permutations (default 100)
        %
        % OUTPUTS:
        %   p_value   - fraction of null runs with accuracy >= observed
        %   null_accs - (n_permutations x 1) null accuracy distribution
        %   obs_acc   - observed mean accuracy from unpermuted classify()
        %
        % EXAMPLE:
        %   [p, null, obs] = Classifier.permutationTest(X, Y, opts, 500);
        %   fprintf('Observed: %.3f, p = %.4f\n', obs, p);
            arguments
                X              table
                Y
                opts           struct = struct()
                n_permutations (1,1) double = 100
            end

            % Observed accuracy
            obs_result = Classifier.classify(X, Y, opts);
            obs_summary = ClassificationResult.summarizeFolds(obs_result);
            obs_acc = obs_summary.mean_accuracy;

            % Permutation null distribution — suppress imbalance warning during shuffles.
            % Each permutation uses a unique derived seed so both Y-shuffle and CV
            % partition vary across runs (avoids conditioning the null on one partition).
            null_accs = zeros(n_permutations, 1);
            perm_opts = opts;
            base_seed = opts.Seed;
            w = warning('off', 'Classifier:imbalanced');
            try
                for p = 1:n_permutations
                    if ~isempty(base_seed)
                        perm_opts.Seed = base_seed + p;
                    end
                    Y_perm = Y(randperm(numel(Y)));
                    perm_result = Classifier.classify(X, Y_perm, perm_opts);
                    perm_summary = ClassificationResult.summarizeFolds(perm_result);
                    null_accs(p) = perm_summary.mean_accuracy;
                end
            finally
                warning(w);
            end

            p_value = mean(null_accs >= obs_acc);

            AnalysisLog.instance().add('Classifier.permutationTest', ...
                struct('n_permutations', n_permutations, 'opts', opts), ...
                sprintf('obs=%.3f p=%.4f (n=%d)', obs_acc, p_value, n_permutations), 0);
        end

        function results = classifyByGroups(X, Y, feature_group_defs, opts)
        % CLASSIFYBYGROUPS  Run classify() over a set of feature-column subsets.
        %
        % feature_group_defs is a struct array where each element has:
        %   .Name    - string label for the group
        %   .Columns - string array of exact column names from X.Properties.VariableNames
        %              (use FeatureStore.unitMatrix(feature_groups) for group-based selection)
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
                sel = ismember(all_cols, cols);
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
        % PREDICT  Apply a trained MATLAB classifier to new data.
        %
        %   trained_model - raw MATLAB classifier (from Classifier.classify result.Mdl)
        %   X             - (N x F) table of feature values
            [Y_pred, scores] = predict(trained_model, X.Variables);
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Static, Access = private)

        function opts = defaultOpts(opts)
            if ~isfield(opts, 'Algorithm'),          opts.Algorithm = 'rf';       end
            if ~isfield(opts, 'KFold'),              opts.KFold = 5;              end
            if ~isfield(opts, 'NHyper'),             opts.NHyper = 0;             end
            if ~isfield(opts, 'CVGroups'),           opts.CVGroups = [];          end
            if ~isfield(opts, 'CVLevel'),            opts.CVLevel = 'recording';  end
            if ~isfield(opts, 'NormalizationGroups'), opts.NormalizationGroups = []; end
            if ~isfield(opts, 'NormalizationVar'),   opts.NormalizationVar = '';  end
            if ~isfield(opts, 'PoolingValues'),      opts.PoolingValues = {};     end
            if ~isfield(opts, 'Prior'),              opts.Prior = 'empirical';    end
            if ~isfield(opts, 'Seed'),               opts.Seed = [];              end
        end

        function checkClassBalance(Y)
            MLUtils.checkClassBalance(Y);
        end

        function [Y_out, group_labels] = prepareLabels(Y, pooling_vals)
        % Convert Y to a consistent format and apply any label pooling.
            if ischar(Y) || (iscell(Y) && ischar(Y{1}))
                Y = string(Y);
            end
            if ~isempty(pooling_vals) && iscell(pooling_vals) && ~isempty(pooling_vals{1})
                [group_labels, ~, group_idx] = unique(Y, 'stable');
                [Y_out, group_labels] = MLPipeline.poolMetadataValues( ...
                    double(group_idx), group_labels, string(pooling_vals));
            else
                [group_labels, ~, Y_out] = unique(Y, 'stable');
            end
        end

        function g = resolveGroups(groups_input, N)
            g = MLUtils.resolveGroups(groups_input, N);
        end

        function [X_train, X_test, feat_names] = normalizeFeatures(X, train_idx, test_idx, norm_groups)
            [X_train, X_test, feat_names] = MLUtils.normalizeFeatures(X, train_idx, test_idx, norm_groups);
        end

    end
end
