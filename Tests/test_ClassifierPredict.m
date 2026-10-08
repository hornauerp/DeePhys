classdef test_ClassifierPredict < matlab.unittest.TestCase
% TEST_CLASSIFIERPREDICT  Fold preprocessing is stored and replayed on new data.
%
% Run: runtests('test_ClassifierPredict')

    methods (Test)

        function testReplayMatchesFoldTransform(tc)
            % applyPreprocessing on the test rows must reproduce normalizeFeatures' X_test.
            rng(1);
            [X, ~, g] = test_ClassifierPredict.makeData(120);
            X{5, 'f2'} = NaN;                         % exercises median imputation
            train_idx = (1:120)' <= 90;
            test_idx  = ~train_idx;
            [~, X_test, ~, prep] = MLUtils.normalizeFeatures(X, train_idx, test_idx, g, "", []);
            X_replay = MLUtils.applyPreprocessing(prep, X(test_idx, :), g(test_idx));
            tc.verifyEqual(X_replay, X_test, 'AbsTol', 1e-12);
        end

        function testGroupStatsKeyedByName(tc)
            % Statistics fitted for plate "p2" must be applied to new "p2" rows even when
            % the new data contains only that plate (indices would map it to the first group).
            rng(2);
            X = [randn(30, 2); randn(30, 2) + 10];
            g = [repmat("p1", 30, 1); repmat("p2", 30, 1)];
            np = NormalizationPipeline.groupThenGlobal();
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            [~, np] = np.fit_transform(X, g);
            x_new = np.transform([10 10], "p2");
            tc.verifyLessThan(max(abs(x_new)), 3);
        end

        function testCellstrGroupsAreNotIgnored(tc)
            labels = MLUtils.groupLabels({'a'; 'b'; 'a'}, 3);
            tc.verifyEqual(labels, ["a"; "b"; "a"]);
            tc.verifyEqual(MLUtils.resolveGroups({'a'; 'b'; 'a'}, 3), [1; 2; 1]);
            tc.verifyEqual(MLUtils.groupLabels(datetime(2024, 1, 1) + [0; 0], 2), ...
                string(datetime(2024, 1, 1) + [0; 0]));
        end

        function testPredictFromResultOnShiftedBatch(tc)
            % New data from an unseen batch with a large offset: replaying group
            % normalisation recovers the classes; the raw model on raw features does not.
            rng(3);
            [X, Y, g] = test_ClassifierPredict.makeData(200);
            opts = struct('Algorithm', 'rf', 'KFold', 2, 'Seed', 1, 'NormalizationGroups', g);
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            result = Classifier.classify(X, Y, opts);
            [X_new, Y_new, g_new] = test_ClassifierPredict.makeData(80);
            X_new.Variables = X_new.Variables + 50;  % batch effect
            g_new = "new_" + g_new;                  % three unseen batches
            Y_pred = Classifier.predict(result(1), X_new, g_new);
            tc.verifyClass(Y_pred, class(Y));
            tc.verifyGreaterThan(mean(Y_pred == Y_new), 0.85);
            tc.verifyError(@() Classifier.predict(result(1), X_new), 'MLUtils:normGroupsRequired');
        end

        function testPredictMatchesColumnsByName(tc)
            rng(4);
            [X, Y] = test_ClassifierPredict.makeData(150);
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            result = Classifier.classify(X, Y, struct('Algorithm', 'rf', 'KFold', 2, 'Seed', 2));
            [X_new, ~] = test_ClassifierPredict.makeData(40);
            p1 = Classifier.predict(result(1), X_new);
            p2 = Classifier.predict(result(1), X_new(:, end:-1:1));
            tc.verifyEqual(p1, p2);
            tc.verifyError(@() Classifier.predict(result(1), X_new(:, 2:end)), 'MLUtils:missingFeatures');
        end

        function testTrainedModelCarriesPreprocessing(tc)
            rng(5);
            [X, Y] = test_ClassifierPredict.makeData(100);
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            result = Classifier.classify(X, Y, struct('Algorithm', 'rf', 'KFold', 2, 'Seed', 3));
            tm = TrainedModel.fromClassificationResult(result(1));
            tc.verifyEqual(tm.FeatureSpec.feature_names, result(1).Preprocessing.InputFeatures);
            tc.verifyNotEmpty(tm.Preprocessing);
        end

        function testTrainedModelPredictsFromTableAndFeatureStore(tc)
            % New-pipeline inputs: a feature table or a FeatureStore (no legacy objects).
            rng(6);
            [X, Y] = test_ClassifierPredict.makeData(120);
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            result = Classifier.classify(X, Y, struct('Algorithm', 'rf', 'KFold', 2, 'Seed', 4));
            tm = TrainedModel.fromClassificationResult(result(1));
            tm.FeatureSpec.level = "Unit";
            [X_new, Y_new] = test_ClassifierPredict.makeData(40);
            p_tab = tm.predict(X_new);
            tc.verifyGreaterThan(mean(p_tab == Y_new), 0.85);
            fs = FeatureStore();
            fs.UnitTable = [table(compose("u%d", (1:40)'), repmat("r1", 40, 1), ...
                'VariableNames', {'UnitID', 'RecordingID'}), X_new];
            fs.MetadataTable = table("r1", 'VariableNames', {'RecordingID'});
            tc.verifyEqual(tm.predict(fs), p_tab);
        end

    end

    methods (Static)

        function [X, Y, g] = makeData(n)
            % Two classes separated on f1/f2; per-batch offsets on all features.
            Y = repmat(["E"; "I"], ceil(n / 2), 1);
            Y = Y(1:n);
            is_i = Y == "I";
            g = repmat(["b1"; "b2"; "b3"], ceil(n / 3), 1);
            g = g(1:n);
            offset = 3 * double(categorical(g));
            f1 = randn(n, 1) + 3 * is_i + offset;
            f2 = randn(n, 1) - 3 * is_i + offset;
            f3 = randn(n, 1) + offset;
            X = table(f1, f2, f3);
        end

    end
end
