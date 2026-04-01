classdef Regressor
% REGRESSOR  Static regression methods that operate on plain feature tables.
%
% Zero coupling to MEArecording, Unit, or RecordingGroup.
%
% USAGE:
%   opts.KFold             = 5;
%   opts.NHyper            = 0;
%   opts.NormalizationGroups = fs.RecordingTable.PlatingDate;
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
        %     .KFold             - number of folds (default 5)
        %     .NHyper            - hyperparam evaluations (default 0)
        %     .NormalizationGroups - (N x 1) string/numeric for per-group z-score
        %     .StratificationVar - (N x 1) for stratified CV splitting
        %
        % OUTPUTS:
        %   result - (1 x K) RegressionResult array
            arguments
                X    table
                Y    {isnumeric}
                opts struct = struct()
            end

            opts = Regressor.defaultOpts(opts);
            norm_groups = Regressor.resolveGroups(opts.NormalizationGroups, height(X));

            % CV partition
            K = opts.KFold;
            if ~isempty(opts.StratificationVar)
                strat = opts.StratificationVar;
                cv = cvpartition(strat, 'KFold', K);
            else
                cv = cvpartition(height(X), 'KFold', K);
            end

            t_start = tic;
            t = templateTree('Surrogate', 'on', 'MinLeafSize', 5, 'NumVariablesToSample', 'auto');
            result(K) = RegressionResult();

            for k = 1:K
                train_idx = cv.training(k);
                test_idx  = cv.test(k);

                Y_train = Y(train_idx);
                Y_test  = Y(test_idx);

                [X_train, X_test] = Regressor.normalizeFeatures(X, train_idx, test_idx, norm_groups);

                clf = fitrensemble(X_train, Y_train, 'Method', 'Bag', ...
                    'NumLearningCycles', 500, 'Learners', t, ...
                    'Options', statset('UseParallel', true));
                Y_pred = predict(clf, X_test);

                result(k) = RegressionResult( ...
                    'Mdl',       clf, ...
                    'Y_pred',    Y_pred, ...
                    'Y_test',    Y_test, ...
                    'mse_train', resubLoss(clf), ...
                    'objects',   [], ...
                    'predImp',   [], ...
                    'Features',  [], ...
                    'Parameters', struct('normalization_var', opts.NormalizationGroups, 'K_fold', K));
            end

            elapsed = toc(t_start);
            AnalysisLog.instance().add('Regressor.regress', opts, ...
                sprintf('%d folds, RMSE=%.4f', K, sqrt(mean([result.mse_train]))), elapsed);
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Static, Access = private)

        function opts = defaultOpts(opts)
            if ~isfield(opts, 'KFold'),              opts.KFold = 5;              end
            if ~isfield(opts, 'NHyper'),             opts.NHyper = 0;             end
            if ~isfield(opts, 'NormalizationGroups'), opts.NormalizationGroups = []; end
            if ~isfield(opts, 'StratificationVar'),  opts.StratificationVar = []; end
        end

        function g = resolveGroups(groups_input, N)
            if isempty(groups_input)
                g = [];
                return
            end
            if ischar(groups_input) || isstring(groups_input)
                [~, ~, g] = unique(string(groups_input(:)), 'stable');
            else
                [~, ~, g] = unique(groups_input(:), 'stable');
            end
            if numel(g) ~= N
                error('Regressor:resolveGroups', ...
                    'Group vector length (%d) must match number of observations (%d).', numel(g), N);
            end
        end

        function [X_train, X_test] = normalizeFeatures(X, train_idx, test_idx, norm_groups)
            mat = X.Variables;

            % Impute NaN
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
        end

    end
end
