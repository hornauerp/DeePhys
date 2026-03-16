classdef MLPipeline < handle
% MLPIPELINE  Machine learning pipeline for classification, regression, and dimensionality reduction.
%
%   Encapsulates the shared workflow: feature assembly → normalization → train/predict → results.
%   Extracted from RecordingGroup to separate ML logic from data container.
%
%   USAGE:
%     mlp = MLPipeline(recording_group);
%     result = mlp.classify("Recording", "rf", "Mutation");
%     result = mlp.regress("Recording", "Concentration");
%     result = mlp.reduceDim("Unit", "UMAP");
%
%   Default hyperparameters are stored in mlp.DefaultParams and can be overridden.

    properties
        DataSource          % RecordingGroup (or any object providing Recordings, Units, Cultures)
        Normalization       % NormalizationPipeline (optional override)
        DefaultParams       % Struct of default hyperparameters
    end

    methods

        function mlp = MLPipeline(data_source, normalization)
        % Constructor.
        %   mlp = MLPipeline(rg)
        %   mlp = MLPipeline(rg, normalization_pipeline)
            arguments
                data_source = []
                normalization = []
            end
            mlp.DataSource = data_source;
            mlp.Normalization = normalization;
            mlp.DefaultParams = MLPipeline.returnDefaultParams();
        end

    end

    methods (Static)

        function params = returnDefaultParams()
        % RETURNDEFAULTPARAMS  Default ML hyperparameters.
        %   Replaces hardcoded values scattered across RecordingGroup methods.
            params.RF.NumCycles = 500;
            params.RF.MinLeafSize = 1;
            params.RF.NumVariablesToSample = 'all';
            params.RF.Surrogate = 'on';
            params.RF.Reproducible = true;

            params.UMAP.NNeighbors = 100;
            params.UMAP.NNeighborsCulture = 10;
            params.UMAP.MinDist = 1;
            params.UMAP.Spread = 5;
            params.UMAP.SGDTasks = 20;
            params.UMAP.ClusterDetail = 'adaptive';

            params.PCA = struct();
            params.tSNE = struct();
        end

        function [clf, train_acc] = createClassifier(X_train, Y_train, alg, N_hyper, params)
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
            end

            if N_hyper > 0
                opt_opts = struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                    'MaxObjectiveEvaluations', N_hyper, 'ShowPlots', false, 'Verbose', 0);
                switch alg
                    case 'svm'
                        clf = fitcsvm(X_train, Y_train, 'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'cnb'
                        clf = fitcnb(X_train, Y_train, 'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'knn'
                        clf = fitcknn(X_train, Y_train, 'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'rf'
                        hyperparams = {'NumLearningCycles', 'MinLeafSize', 'MaxNumSplits', 'SplitCriterion', 'NumVariablesToSample'};
                        t = templateTree('Reproducible', true);
                        clf = fitcensemble(X_train, Y_train, 'Method', 'Bag', ...
                            'OptimizeHyperparameters', hyperparams, 'Learners', t, ...
                            'HyperparameterOptimizationOptions', opt_opts);
                end
                train_acc = 1 - clf.HyperparameterOptimizationResults.MinObjective;  % MinObjective is error rate → convert to accuracy
            else
                switch alg
                    case 'svm'
                        clf = fitcsvm(X_train, Y_train);
                    case 'cnb'
                        clf = fitcnb(X_train, Y_train);
                    case 'knn'
                        clf = fitcknn(X_train, Y_train);
                    case 'rf'
                        t = templateTree('Surrogate', params.RF.Surrogate, ...
                            'MinLeafSize', params.RF.MinLeafSize, ...
                            'NumVariablesToSample', params.RF.NumVariablesToSample, ...
                            'Reproducible', params.RF.Reproducible);
                        clf = fitcensemble(X_train, Y_train, 'Method', 'Bag', ...
                            'NumLearningCycles', params.RF.NumCycles, ...
                            'Learners', t, 'Options', statset("UseParallel", true));
                end
                train_acc = 1 - resubLoss(clf, 'LossFun', 'classerror');
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
            clf_group_idx = find(contains(group_labels, classification_val));
            new_group_idx = ismember(group_idx, clf_group_idx) * 1;
            new_group_labels(1) = join(group_labels(clf_group_idx), '/');
            clf_group_idx = find(~ismember(group_labels, classification_val));
            new_group_idx(new_group_idx == 0) = 2;
            new_group_labels(2) = join(group_labels(clf_group_idx), '/');
        end

    end
end
