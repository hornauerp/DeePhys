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

                [X_train, X_test, feat_names] = Regressor.normalizeFeatures(X, train_idx, test_idx, norm_groups);

                [mdl, train_r2] = MLPipeline.createRegressor(X_train, Y_train, opts.Algorithm, opts.NHyper);
                Y_pred = predict(mdl, X_test);

                if opts.Algorithm == "rf"
                    predImp = oobPermutedPredictorImportance(mdl, 'Options', statset('UseParallel', true));
                else
                    predImp = [];
                end

                result(k) = RegressionResult( ...
                    'Mdl',       mdl, ...
                    'Y_pred',    Y_pred, ...
                    'Y_test',    Y_test, ...
                    'mse_train', resubLoss(mdl), ...
                    'train_acc', train_r2, ...
                    'objects',   [], ...
                    'predImp',   predImp, ...
                    'Features',  feat_names, ...
                    'Parameters', struct( ...
                        'algorithm',        opts.Algorithm, ...
                        'normalization_var', opts.NormalizationGroups, ...
                        'K_fold',           K, ...
                        'N_hyper',          opts.NHyper, ...
                        'seed',             opts.Seed));
            end

            elapsed = toc(t_start);
            summary = RegressionResult.summarizeFolds(result);
            AnalysisLog.instance().add('Regressor.regress', opts, ...
                sprintf('%d folds, %s, R2=%.3f', K, opts.Algorithm, summary.mean_R2), elapsed);
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
            if ~isfield(opts, 'Seed'),               opts.Seed = [];              end
            if ~isfield(opts, 'NormalizationGroups'), opts.NormalizationGroups = []; end
            if ~isfield(opts, 'CVGroups'),           opts.CVGroups = [];          end
            if ~isfield(opts, 'StratificationVar'),  opts.StratificationVar = []; end
        end

        function g = resolveGroups(groups_input, N)
            g = MLUtils.resolveGroups(groups_input, N);
        end

        function [X_train, X_test, feat_names] = normalizeFeatures(X, train_idx, test_idx, norm_groups)
            mat        = X.Variables;
            feat_names = string(X.Properties.VariableNames);

            % Impute NaN using training-set median
            for col = 1:size(mat, 2)
                train_col = mat(train_idx, col);
                med = median(train_col(~isnan(train_col)));
                if isnan(med), med = 0; end
                mat(isnan(mat(:, col)), col) = med;
            end

            if ~isempty(norm_groups)
                np = NormalizationPipeline.groupThenGlobal();
                g_train = norm_groups(train_idx);
            else
                np = NormalizationPipeline.globalOnly();
                g_train = [];
            end

            [X_train, np] = np.fit_transform(mat(train_idx, :), g_train);

            if ~isempty(norm_groups)
                X_test = np.transform(mat(test_idx, :), norm_groups(test_idx));
            else
                X_test = np.transform(mat(test_idx, :), []);
            end

            bad = any(isnan(X_train)) | any(isnan(X_test));
            X_train(:, bad) = [];
            X_test(:, bad)  = [];
            feat_names(bad) = [];
        end

    end
end
