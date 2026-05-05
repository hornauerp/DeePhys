classdef Regressor
% REGRESSOR  Static regression methods that operate on plain feature tables.
%
% Zero coupling to MEArecording, Unit, or RecordingGroup.
%
% USAGE:
%   opts.Algorithm         = 'rf';          % 'rf', 'svm', 'knn'
%   opts.KFold             = 5;
%   opts.NHyper            = 0;             % 0 = no hyperparam search
%   opts.Seed              = 42;            % RNG seed for reproducibility
%   opts.NormalizationGroups = fs.RecordingTable.PlatingDate;
%   opts.CVGroups            = fs.UnitTable.RecordingID;  % hierarchy-aware CV
%
%   result = Regressor.regress(X, Y, opts);
%   [p, null, obs] = Regressor.permutationTest(X, Y, opts, 500);
%   results = Regressor.regressByGroups(X, Y, feature_group_defs, opts);

    methods (Static)

        function result = regress(X, Y, opts)
        % REGRESS  K-fold cross-validated regression on a feature table.
        %
        % INPUTS:
        %   X    - (N x F) table of feature values
        %   Y    - (N x 1) numeric target vector
        %   opts - struct with optional fields:
        %     .Algorithm         - 'rf' (default), 'svm', 'knn'
        %     .KFold             - number of folds (default 5)
        %     .NHyper            - Bayesian hyperparam evaluations (default 0)
        %     .Seed              - RNG seed for reproducibility (default [])
        %     .NormalizationGroups - (N x 1) string/numeric for per-group z-score
        %     .NormalizationVar    - metadata field name for serialization (default '')
        %     .CVGroups          - (N x 1) string/numeric for hierarchy-aware CV
        %                          (e.g. fs.UnitTable.RecordingID). [] = standard CV.
        %     .StratificationVar - (N x 1) for stratified CV splitting (ignored when CVGroups set)
        %
        % OUTPUTS:
        %   result - (1 x K) RegressionResult array
            arguments
                X    table
                Y    {isnumeric}
                opts struct = struct()
            end

            opts = Regressor.defaultOpts(opts);

            % Reproducibility
            if ~isempty(opts.Seed)
                rng(opts.Seed);
            end

            norm_groups = Regressor.resolveGroups(opts.NormalizationGroups, height(X));
            cv_groups   = Regressor.resolveGroups(opts.CVGroups, height(X));

            % CV partition — prefer hierarchy-aware GroupedCV when CVGroups provided
            K = opts.KFold;
            if ~isempty(cv_groups)
                gcv     = GroupedCV.byGroups(cv_groups, K);
                K       = gcv.NumFolds;
                use_gcv = true;
            elseif ~isempty(opts.StratificationVar)
                cv      = cvpartition(opts.StratificationVar, 'KFold', K);
                use_gcv = false;
            else
                cv      = cvpartition(height(X), 'KFold', K);
                use_gcv = false;
            end

            t_start = tic;
            valid_folds = true(1, K);
            result(K) = RegressionResult();

            for k = 1:K
                if use_gcv
                    [train_idx, test_idx] = gcv.fold(k);
                else
                    train_idx = cv.training(k);
                    test_idx  = cv.test(k);
                end

                Y_train = Y(train_idx);
                Y_test  = Y(test_idx);

                [X_train, X_test, feat_names] = Regressor.normalizeFeatures(X, train_idx, test_idx, norm_groups, opts.NormalizationPipeline, opts.CovariateLabels);

                if isempty(X_test)
                    warning('Regressor:emptyTestFold', ...
                        'Fold %d has an empty test set — skipped.', k);
                    valid_folds(k) = false;
                    continue
                end

                [mdl, train_r2] = MLPipeline.createRegressor(X_train, Y_train, opts.Algorithm, opts.NHyper);
                Y_pred = predict(mdl, X_test);

                if opts.Algorithm == "rf"
                    predImp = oobPermutedPredictorImportance(mdl, 'Options', statset('UseParallel', true));
                else
                    predImp = [];
                end

                if opts.Algorithm == "rf"
                    mse_train_val = oobLoss(mdl);       % OOB ≈ held-out performance, not true training MSE
                else
                    mse_train_val = resubLoss(mdl);     % resubstitution; optimistic (~0 for SVM/KNN)
                end
                result(k) = RegressionResult( ...
                    'Mdl',       mdl, ...
                    'Y_pred',    Y_pred, ...
                    'Y_test',    Y_test, ...
                    'mse_train', mse_train_val, ...
                    'train_acc', train_r2, ...
                    'objects',   [], ...
                    'predImp',   predImp, ...
                    'Features',  feat_names, ...
                    'Parameters', struct( ...
                        'algorithm',        opts.Algorithm, ...
                        'normalization_var', opts.NormalizationVar, ...
                        'K_fold',           K, ...
                        'N_hyper',          opts.NHyper, ...
                        'seed',             opts.Seed));
            end
            result = result(valid_folds);

            elapsed = toc(t_start);
            summary = RegressionResult.summarizeFolds(result);
            AnalysisLog.instance().add('Regressor.regress', opts, ...
                sprintf('%d folds, %s, R2=%.3f', K, opts.Algorithm, summary.mean_R2), elapsed);
        end

        function [p_value, null_r2s, obs_r2] = permutationTest(X, Y, opts, n_permutations)
        % PERMUTATIONTEST  Estimate regression significance via label permutation.
        %
        % Shuffles target values n_permutations times and runs regress() each time,
        % building a null distribution of R² values. Returns a one-sided p-value:
        %   P(null_R² >= observed_R²)
        %
        % INPUTS:
        %   X              - (N x F) feature table
        %   Y              - (N x 1) numeric target vector
        %   opts           - same opts struct as Regressor.regress (reused each run)
        %   n_permutations - number of permutations (default 100)
        %
        % OUTPUTS:
        %   p_value   - fraction of null runs with R² >= observed
        %   null_r2s  - (n_permutations x 1) null R² distribution
        %   obs_r2    - observed mean R² from unpermuted regress()
            arguments
                X              table
                Y              {isnumeric}
                opts           struct = struct()
                n_permutations (1,1) double = 100
            end

            obs_result  = Regressor.regress(X, Y, opts);
            obs_summary = RegressionResult.summarizeFolds(obs_result);
            obs_r2      = obs_summary.mean_R2;

            null_r2s  = zeros(n_permutations, 1);
            perm_opts = opts;
            base_seed = opts.Seed;
            for p = 1:n_permutations
                if ~isempty(base_seed)
                    perm_opts.Seed = base_seed + p;
                    rng(base_seed + p);  % seed the Y-shuffle too, not just the CV partition inside regress()
                end
                Y_perm = Y(randperm(numel(Y)));
                perm_result  = Regressor.regress(X, Y_perm, perm_opts);
                perm_summary = RegressionResult.summarizeFolds(perm_result);
                null_r2s(p)  = perm_summary.mean_R2;
            end

            p_value = mean(null_r2s >= obs_r2);

            AnalysisLog.instance().add('Regressor.permutationTest', ...
                struct('n_permutations', n_permutations, 'opts', opts), ...
                sprintf('obs_R2=%.3f p=%.4f (n=%d)', obs_r2, p_value, n_permutations), 0);
        end

        function results = regressByGroups(X, Y, feature_group_defs, opts)
        % REGRESSBYGROUPS  Run regress() over a set of feature-column subsets.
        %
        % feature_group_defs is a struct array where each element has:
        %   .Name    - string label for the group
        %   .Columns - string array of exact column names from X.Properties.VariableNames
        %
        % Returns a struct array with .Name and .Result fields.
            arguments
                X                  table
                Y                  {isnumeric}
                feature_group_defs struct
                opts               struct = struct()
            end
            results = struct('Name', {}, 'Result', {});
            for g = 1:numel(feature_group_defs)
                fg   = feature_group_defs(g);
                cols = string(fg.Columns);
                all_cols = string(X.Properties.VariableNames);
                sel  = ismember(all_cols, cols);
                if ~any(sel)
                    warning('Regressor:regressByGroups', ...
                        'No columns matched for group "%s" — skipped.', fg.Name);
                    continue
                end
                Xsub = X(:, sel);
                res  = Regressor.regress(Xsub, Y, opts);
                results(end+1) = struct('Name', fg.Name, 'Result', res); %#ok<AGROW>
            end
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Static, Access = private)

        function opts = defaultOpts(opts)
            if ~isfield(opts, 'Algorithm'),           opts.Algorithm = 'rf';       end
            if ~isfield(opts, 'KFold'),               opts.KFold = 5;              end
            if ~isfield(opts, 'NHyper'),              opts.NHyper = 0;             end
            if ~isfield(opts, 'Seed'),                opts.Seed = [];              end
            if ~isfield(opts, 'NormalizationGroups'),  opts.NormalizationGroups = [];  end
            if ~isfield(opts, 'NormalizationVar'),    opts.NormalizationVar = '';     end
            if ~isfield(opts, 'NormalizationPipeline'), opts.NormalizationPipeline = ''; end
            if ~isfield(opts, 'CovariateLabels'),     opts.CovariateLabels = [];     end
            if ~isfield(opts, 'CVGroups'),            opts.CVGroups = [];          end
            if ~isfield(opts, 'StratificationVar'),   opts.StratificationVar = []; end
        end

        function g = resolveGroups(groups_input, N)
            g = MLUtils.resolveGroups(groups_input, N);
        end

        function [X_train, X_test, feat_names] = normalizeFeatures(X, train_idx, test_idx, norm_groups, pipeline_type, covariate_labels)
            if nargin < 5, pipeline_type = ""; end
            if nargin < 6, covariate_labels = []; end
            [X_train, X_test, feat_names] = MLUtils.normalizeFeatures(X, train_idx, test_idx, norm_groups, pipeline_type, covariate_labels);
        end

    end
end
