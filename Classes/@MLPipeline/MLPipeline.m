classdef MLPipeline
% MLPIPELINE  Static ML helpers: classifier/regressor creation, CV splits, label pooling.
%
%   All methods are static — no instantiation needed.
%   Access default hyperparameters via MLPipeline.returnDefaultParams().

    methods (Static)

        function params = returnDefaultParams()
        % RETURNDEFAULTPARAMS  Default ML hyperparameters.
        %   Replaces hardcoded values scattered across RecordingGroup methods.
            params.RF.NumCycles = 500;
            params.RF.MinLeafSize = 5;  % ≥5 prevents leaf-level memorisation on typical neuro datasets
            params.RF.NumVariablesToSample = 'auto';  % auto = sqrt(F) for classif, F/3 for regress (Breiman, 2001)
            params.RF.Surrogate = 'on';
            params.RF.Reproducible = true;
            params.RF.Prior = 'empirical';  % 'uniform' gives equal weight to each class regardless of frequency

            params.UMAP.NNeighbors = 100;
            params.UMAP.NNeighborsCulture = 10;
            params.UMAP.MinDist = 1;
            params.UMAP.Spread = 5;
            params.UMAP.SGDTasks = 20;
            params.UMAP.ClusterDetail = 'adaptive';

            params.PCA = struct();
            params.tSNE = struct();
        end

        function [clf, train_acc] = createClassifier(X_train, Y_train, alg, N_hyper, params, train_groups)
        % CREATECLASSIFIER  Train a classifier with optional hyperparameter optimization.
        %
        % INPUTS:
        %   X_train - (N x F) training feature matrix
        %   Y_train - (N x 1) training labels
        %   alg     - "rf", "svm", "cnb", "knn"
        %   N_hyper - number of hyperparameter optimization evaluations (0 = none)
        %   params  - (optional) struct with RF hyperparameters
            arguments
                X_train {isnumeric}
                Y_train
                alg string = "rf"
                N_hyper (1,1) double = 0
                params struct = MLPipeline.returnDefaultParams()
                train_groups = []   % CV groups of the training rows: grouped inner CV
            end

            if N_hyper > 0
                opt_opts = struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                    'MaxObjectiveEvaluations', N_hyper, 'ShowPlots', false, 'Verbose', 0);
                % Inner folds keep each group (culture) together; random folds reward
                % hyperparameters that memorise groups
                if ~isempty(train_groups)
                    cvp = MLUtils.groupedPartition(train_groups, 5);
                    if ~isempty(cvp), opt_opts.CVPartition = cvp; end
                end
                switch alg
                    case 'svm'
                        if numel(unique(Y_train)) > 2   % fitcsvm is binary-only
                            clf = fitcecoc(X_train, Y_train, 'Learners', templateSVM(), 'Prior', params.RF.Prior, ...
                                'OptimizeHyperparameters', 'auto', ...
                                'HyperparameterOptimizationOptions', opt_opts);
                        else
                            clf = fitcsvm(X_train, Y_train, 'Prior', params.RF.Prior, ...
                                'OptimizeHyperparameters', 'all', ...
                                'HyperparameterOptimizationOptions', opt_opts);
                        end
                    case 'cnb'
                        clf = fitcnb(X_train, Y_train, 'Prior', params.RF.Prior, ...
                            'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'knn'
                        clf = fitcknn(X_train, Y_train, 'Prior', params.RF.Prior, ...
                            'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'rf'
                        hyperparams = {'NumLearningCycles', 'MinLeafSize', 'MaxNumSplits', 'SplitCriterion', 'NumVariablesToSample'};
                        t = templateTree('Reproducible', true);
                        clf = fitcensemble(X_train, Y_train, 'Method', 'Bag', 'Prior', params.RF.Prior, ...
                            'OptimizeHyperparameters', hyperparams, 'Learners', t, ...
                            'HyperparameterOptimizationOptions', opt_opts);
                end
                train_acc = 1 - clf.HyperparameterOptimizationResults.MinObjective;  % MinObjective is error rate → convert to accuracy
            else
                switch alg
                    case 'svm'
                        if numel(unique(Y_train)) > 2
                            % fitcsvm is binary-only; one-vs-one ECOC for multiclass labels
                            clf = fitcecoc(X_train, Y_train, 'Learners', templateSVM(), 'Prior', params.RF.Prior);
                        else
                            clf = fitcsvm(X_train, Y_train, 'Prior', params.RF.Prior);
                        end
                    case 'cnb'
                        clf = fitcnb(X_train, Y_train, 'Prior', params.RF.Prior);
                    case 'knn'
                        clf = fitcknn(X_train, Y_train, 'Prior', params.RF.Prior);
                    case 'rf'
                        t = templateTree('Surrogate', params.RF.Surrogate, ...
                            'MinLeafSize', params.RF.MinLeafSize, ...
                            'NumVariablesToSample', MLPipeline.resolveNumVariables( ...
                                params.RF.NumVariablesToSample, size(X_train, 2), "classification"), ...
                            'Reproducible', params.RF.Reproducible);
                        clf = fitcensemble(X_train, Y_train, 'Method', 'Bag', 'Prior', params.RF.Prior, ...
                            'NumLearningCycles', params.RF.NumCycles, ...
                            'Learners', t, 'Options', statset("UseParallel", true));
                end
                if alg == "rf"
                    % OOB accuracy: held-out only for exchangeable rows; with grouped rows (units
                    % of a culture) the OOB rows share groups with the in-bag rows -> optimistic
                    train_acc = 1 - oobLoss(clf, 'LossFun', 'classiferror');
                else
                    train_acc = 1 - resubLoss(clf, 'LossFun', 'classiferror');  % resubstitution; optimistic (~1.0 for SVM/KNN)
                end
            end
        end

        function n = resolveNumVariables(spec, n_features, task)
        % RESOLVENUMVARIABLES  Predictors sampled per split. "auto" follows Breiman (2001):
        %   sqrt(F) for classification, F/3 for regression. templateTree itself only
        %   accepts a positive integer or 'all' ('all' = plain bagging, not a random forest).
            if (ischar(spec) || isstring(spec)) && string(spec) == "auto"
                if task == "classification"
                    n = max(1, round(sqrt(n_features)));
                else
                    n = max(1, round(n_features / 3));
                end
            else
                n = spec;
            end
        end

        function [mdl, train_r2] = createRegressor(X_train, Y_train, alg, N_hyper, train_groups)
        % CREATEREGRESSOR  Train a regression model with optional hyperparameter optimization.
        %
        % INPUTS:
        %   X_train - (N x F) training feature matrix
        %   Y_train - (N x 1) numeric training targets
        %   alg     - "rf", "svm", "knn"
        %   N_hyper - Bayesian optimization evaluations (0 = none)
        %
        % OUTPUTS:
        %   mdl      - trained MATLAB regression model
        %   train_r2 - training-set R² (in-bag for RF, resubstitution for SVM/KNN)
            arguments
                X_train {isnumeric}
                Y_train {isnumeric}
                alg     string = "rf"
                N_hyper (1,1) double = 0
                train_groups = []   % CV groups of the training rows: grouped inner CV
            end

            params = MLPipeline.returnDefaultParams();
            opt_opts = struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                'MaxObjectiveEvaluations', N_hyper, 'ShowPlots', false, 'Verbose', 0);
            if N_hyper > 0 && ~isempty(train_groups)
                cvp = MLUtils.groupedPartition(train_groups, 5);   % inner folds keep groups together
                if ~isempty(cvp), opt_opts.CVPartition = cvp; end
            end

            if N_hyper > 0
                switch alg
                    case 'svm'
                        mdl = fitrsvm(X_train, Y_train, ...
                            'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'knn'   % MATLAB has no fitrknn: K chosen by leave-one-(group-)out error
                        mdl = KNNRegressionModel(X_train, Y_train, ...
                            KNNRegressionModel.selectK(X_train, Y_train, [], train_groups));
                    case 'rf'
                        hyperparams = {'NumLearningCycles', 'MinLeafSize', 'MaxNumSplits', 'NumVariablesToSample'};
                        t = templateTree('Reproducible', true);
                        mdl = fitrensemble(X_train, Y_train, 'Method', 'Bag', ...
                            'Learners', t, ...
                            'OptimizeHyperparameters', hyperparams, ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    otherwise
                        error('MLPipeline:createRegressor', 'Unknown algorithm "%s". Use rf, svm, or knn.', alg);
                end
            else
                switch alg
                    case 'svm'
                        mdl = fitrsvm(X_train, Y_train);
                    case 'knn'   % MATLAB has no fitrknn
                        mdl = KNNRegressionModel(X_train, Y_train, 5);
                    case 'rf'
                        t = templateTree('Surrogate', params.RF.Surrogate, ...
                            'MinLeafSize', params.RF.MinLeafSize, ...
                            'NumVariablesToSample', MLPipeline.resolveNumVariables( ...
                                params.RF.NumVariablesToSample, size(X_train, 2), "regression"), ...
                            'Reproducible', params.RF.Reproducible);
                        mdl = fitrensemble(X_train, Y_train, 'Method', 'Bag', ...
                            'NumLearningCycles', params.RF.NumCycles, ...
                            'Learners', t, 'Options', statset('UseParallel', true));
                    otherwise
                        error('MLPipeline:createRegressor', 'Unknown algorithm "%s". Use rf, svm, or knn.', alg);
                end
            end

            % Training R² (OOB estimate for RF; resubstitution for SVM/KNN)
            if alg == "rf"
                Y_train_pred = oobPredict(mdl);
            else
                Y_train_pred = predict(mdl, X_train);
            end
            ss_res = sum((Y_train - Y_train_pred).^2);
            ss_tot = sum((Y_train - mean(Y_train)).^2);
            if ss_tot == 0
                train_r2 = NaN;
            else
                train_r2 = 1 - ss_res / ss_tot;
            end
        end

        function [Y_train, Y_test, train_idx, test_idx] = cvSplit(Y, cv, k)
        % CVSPLIT  Extract train/test split for fold k.
            arguments
                Y
                cv cvpartition
                k double
            end
            if iscell(Y)
                Y_train = [Y{cv.training(k)}];
                Y_test  = [Y{cv.test(k)}];
                k_train = cv.training(k);
                train_idx = arrayfun(@(x) ones(size(Y{x})) * k_train(x), 1:length(Y), 'un', 0);
                train_idx = [train_idx{:}];
                test_idx = ~train_idx;
            else
                Y_train = Y(cv.training(k));
                Y_test  = Y(cv.test(k));
                train_idx = cv.training(k);
                test_idx  = cv.test(k);
            end
        end

        function [new_group_idx, new_group_labels] = poolMetadataValues(group_idx, group_labels, classification_val)
        % POOLMETADATAVALUES  Pool metadata values into binary groups.
            arguments
                group_idx double
                group_labels string
                classification_val string
            end
            % Exact matches: a substring match would also pool e.g. "HET_KD" with "HET"
            in_pool = ismember(group_labels, classification_val);
            new_group_idx = ismember(group_idx, find(in_pool)) * 1;
            new_group_labels(1) = join(group_labels(in_pool), '/');
            new_group_idx(new_group_idx == 0) = 2;
            new_group_labels(2) = join(group_labels(~in_pool), '/');
        end

    end
end
