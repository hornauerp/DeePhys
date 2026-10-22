classdef test_ClassifierRFOptions < matlab.unittest.TestCase
% TEST_CLASSIFIERRFOPTIONS  opts.RF overrides of the random forest tree size.
%
% Run: runtests('test_ClassifierRFOptions')

    methods (Test)

        function testTemplateTreeSize(tc)
            % Default template (MinParentSize unset -> templateTree's 10): trees on 8 rows
            % stay a single leaf. MinLeafSize 1 / MinParentSize 2 reach the fitted trees.
            rf = MLPipeline.returnDefaultParams().RF;
            tc.verifyEmpty(rf.MinParentSize);
            [X, Y] = test_ClassifierRFOptions.makeData();
            rng(1);
            mdl = fitcensemble(X{:, :}, Y, 'Method', 'Bag', 'NumLearningCycles', 10, ...
                'Learners', MLPipeline.treeTemplate(rf, 3, "classification"));
            tc.verifyTrue(all(cellfun(@(t) t.NumNodes, mdl.Trained) == 1));
            rf.MinLeafSize = 1; rf.MinParentSize = 2;
            mdl = fitcensemble(X{:, :}, Y, 'Method', 'Bag', 'NumLearningCycles', 10, ...
                'Learners', MLPipeline.treeTemplate(rf, 3, "classification"));
            tc.verifyTrue(any(cellfun(@(t) t.NumNodes, mdl.Trained) > 1));
        end

        function testSmallTreesLearnFewRowsPerClass(tc)
            % Three training rows per class: the default trees cannot split (MinParentSize 10),
            % so every fold predicts one class; MinLeafSize 1 / MinParentSize 2 separates them.
            [X, Y] = test_ClassifierRFOptions.makeData();
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            base = struct('Algorithm', 'rf', 'KFold', 4, 'Seed', 1, 'Prior', 'uniform', ...
                'ComputeImportance', false);
            r_def = Classifier.classify(X, Y, base);
            for k = 1:numel(r_def)
                tc.verifyNumElements(unique(r_def(k).Y_pred), 1);
            end
            tc.verifyEqual(ClassificationResult.summarizeFolds(r_def).balanced_accuracy, 0.5, 'AbsTol', 1e-12);

            opts = base;
            opts.RF = struct('MinLeafSize', 1, 'MinParentSize', 2);
            r_small = Classifier.classify(X, Y, opts);
            tc.verifyEqual(ClassificationResult.summarizeFolds(r_small).balanced_accuracy, 1);
        end

        function testInvalidRFOptions(tc)
            [X, Y] = test_ClassifierRFOptions.makeData();
            base = struct('Algorithm', 'rf', 'KFold', 4, 'ComputeImportance', false);
            bad = {struct('MinLeaf', 1), ...                % typo: unknown field
                   struct('Prior', 'uniform'), ...          % prior belongs in opts.Prior
                   struct('MinParentSize', 0), ...
                   struct('MinLeafSize', 1.5)};
            for i = 1:numel(bad)
                opts = base; opts.RF = bad{i};
                tc.verifyError(@() Classifier.classify(X, Y, opts), 'Classifier:rfOptions');
            end
        end

    end

    methods (Static)

        function [X, Y] = makeData()
            % 2 classes x 4 rows; f1 separates the classes, f2/f3 are noise
            rng(7);
            Y  = categorical(repelem(["A"; "B"], 4));
            f1 = double(Y == "B") * 10 + 0.1 * randn(8, 1);
            X  = table(f1, randn(8, 1), randn(8, 1), 'VariableNames', ["f1", "f2", "f3"]);
        end

    end
end
