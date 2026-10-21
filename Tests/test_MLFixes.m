classdef test_MLFixes < matlab.unittest.TestCase
% TEST_MLFIXES  Leakage, metric and statistics checks of the ML layer on synthetic data.
%
% Features of the "culture fingerprint" data identify the culture only; labels / targets are
% set per culture and unrelated to the features, so any score above chance is leakage.
%
% Run: runtests('test_MLFixes')

    methods (Test)

        function testExperimentRegressGroupsByCulture(tc)
            [ex, ~] = test_MLFixes.fingerprintExperiment(16, 4);
            r = ex.regress("Recording", "Conc", struct('KFold', 5, 'Seed', 1, 'ComputeImportance', false));
            s = RegressionResult.summarizeFolds(r);
            tc.verifyLessThan(s.R2_pooled, 0.1);   % row-level CV gave about +0.3 on this data
        end

        function testPooledR2DefinedForCultureConstantTarget(tc)
            rng(2); [X, cul] = test_MLFixes.fingerprintMatrix(16, 4);
            Y = repelem([1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4], 4);
            r = Regressor.regress(X, Y, struct('KFold', -1, 'CVGroups', cul, 'Seed', 1, 'ComputeImportance', false));
            s = RegressionResult.summarizeFolds(r);
            tc.verifyTrue(all(isnan(s.per_fold.R2)));   % each fold holds one culture, one target value
            tc.verifyFalse(isnan(s.R2_pooled));
        end

        function testPermutationTestRejectsUndefinedR2(tc)
            rng(3); [X, cul] = test_MLFixes.fingerprintMatrix(6, 3);
            Y = ones(height(X), 1);
            tc.verifyError(@() Regressor.permutationTest(X, Y, ...
                struct('KFold', 3, 'CVGroups', cul, 'Seed', 1), 3), 'Regressor:undefinedR2');
        end

        function testRegressDropsNaNTargets(tc)
            rng(4); [X, cul] = test_MLFixes.fingerprintMatrix(8, 4);
            Y = randn(height(X), 1); Y(3) = NaN;
            w = warning('off', 'MLUtils:missingTargets'); c = onCleanup(@() warning(w));
            r = Regressor.regress(X, Y, struct('KFold', 4, 'CVGroups', cul, 'Seed', 1, 'ComputeImportance', false));
            tc.verifyEqual(sum(arrayfun(@(z) numel(z.Y_test), r)), height(X) - 1);
        end

        function testMissingLabelsAreNotAClass(tc)
            rng(5); [X, cul] = test_MLFixes.fingerprintMatrix(12, 10);
            Y = repelem(repmat(["WT"; "KO"], 6, 1), 10); Y(1:3) = missing;
            w = warning('off', 'MLUtils:missingTargets'); c = onCleanup(@() warning(w));
            r = Classifier.classify(X, Y, struct('CVGroups', cul, 'Seed', 1, 'ComputeImportance', false));
            tc.verifyEqual(numel(r(1).GroupLabels), 2);
        end

        function testScoreColumnsFollowClasses(tc)
            rng(6); [X, cul] = test_MLFixes.fingerprintMatrix(12, 10);
            Y = repelem([1;1;1;1;1;2;2;2;2;2;2;3], 10);   % class 3: one culture, absent from a training fold
            r = Classifier.classify(X, Y, struct('CVGroups', cul, 'Seed', 1, 'ComputeImportance', false));
            tc.verifyTrue(all(arrayfun(@(z) size(z.scores, 2), r) == 3));
            T = ClassificationResult.results2table(r);
            tc.verifyEqual(sum(startsWith(string(T.Properties.VariableNames), "Score_")), 3);
        end

        function testFoldF1IsAlwaysMacro(tc)
            r = ClassificationResult('Y_pred', [2;2;2;2], 'Y_test', [2;2;2;3]);
            m = r.computeMetrics();
            % class 2: P 0.75, R 1 -> F1 0.857; class 3: F1 0 -> macro 0.429
            tc.verifyEqual(m.F1_score, (2*0.75/1.75) / 2, 'AbsTol', 1e-12);
        end

        function testTestFoldImportanceIgnoresFingerprints(tc)
            rng(7); nC = 16; nU = 30; F = 10; cul = repelem((1:nC)', nU);
            fp = randn(nC, F); Xf = fp(cul, :) + 0.5 * randn(nC*nU, F);
            Y = repelem(randi(2, nC, 1), nU);
            sig = 1.5 * (Y == 2) + randn(nC*nU, 1);
            X = array2table([Xf, sig], 'VariableNames', ["Fp" + (1:F), "Signal"]);
            r = Classifier.classify(X, Y, struct('CVGroups', cul, 'Seed', 1));
            T = MLResult.summarizeImportance(r);
            tc.verifyEqual(T.Feature(1), "Signal");
        end

        function testKNNSelectKLeaveOneGroupOut(tc)
            rng(8); [X, cul] = test_MLFixes.fingerprintMatrix(10, 20);
            Y = repelem(randn(10, 1), 20);
            K_logo = KNNRegressionModel.selectK(X.Variables, Y, [], cul);
            K_loo  = KNNRegressionModel.selectK(X.Variables, Y);
            tc.verifyLessThanOrEqual(K_loo, K_logo);   % leave-one-out rewards culture mates
            tc.verifyGreaterThan(K_logo, 3);
        end

        function testPairedComparisonMatchesPostIDs(tc)
            pre = (1:8)'; post_ids = (8:-1:1)';
            post = flipud(pre) + 0.5 + 0.01 * (1:8)';
            res = StatisticalTest.pairedComparison(pre, post, [(1:8)'; post_ids], 'test', "ttest");
            tc.verifyEqual(res.mean_diff, 0.5 + 0.01 * mean(1:8), 'AbsTol', 1e-12);
            tc.verifyLessThan(res.p, 1e-6);
        end

        function testCompareGroupsIgnoresNaN(tc)
            rng(9); a = [randn(20, 1); NaN]; b = randn(20, 1) + 1;
            res = StatisticalTest.compareGroups([a; b], [ones(21, 1); 2 * ones(20, 1)]);
            tc.verifyEqual(res.group_ns, [20 20]);
            tc.verifyFalse(isnan(res.effect_size));
            tc.verifyNotEqual(res.effect_size, 0);
        end

        function testMixedEffectsBetweenClusterDF(tc)
            rng(10); nC = 6; nU = 100; c = repelem((1:nC)', nU);
            y = randn(nC, 1); y = y(c) + randn(nC*nU, 1);
            T = table(y, categorical(double(c > 3)), categorical(c), 'VariableNames', {'Feature', 'Group', 'CultureID'});
            R = StatisticalTest.mixedEffects(T, 'Feature ~ Group + (1|CultureID)');
            df = R.anova_table.DF2(2);
            tc.verifyLessThan(df, 10);   % ~ cultures, not ~ units (residual DF would be 598)
        end

    end

    methods (Static)

        function [X, cul] = fingerprintMatrix(nC, nR)
            F = 20; cul = repelem((1:nC)', nR);
            fp = randn(nC, F);
            X = array2table(fp(cul, :) + 0.3 * randn(nC*nR, F));
        end

        function [ex, cul] = fingerprintExperiment(nC, nR)
            rng(1); F = 20; cul = repelem((1:nC)', nR);
            rec = "r" + (1:nC*nR)';
            fp = randn(nC, F); Xf = fp(cul, :) + 0.3 * randn(nC*nR, F);
            yC = randn(nC, 1) * 10;
            MetadataTable = table(rec, "C" + cul, ones(nC*nR, 1), yC(cul), ...
                'VariableNames', {'RecordingID', 'ChipID', 'PlatingDate', 'Conc'});
            RecordingTable = [MetadataTable, array2table(Xf, 'VariableNames', "F" + (1:F))]; %#ok<NASGU>
            UnitTable = table(string.empty(0, 1), string.empty(0, 1), 'VariableNames', {'UnitID', 'RecordingID'}); %#ok<NASGU>
            MetadataFields = ["ChipID", "PlatingDate", "Conc"]; %#ok<NASGU>
            fn = [tempname '.mat'];
            save(fn, 'UnitTable', 'RecordingTable', 'MetadataTable', 'MetadataFields');
            fs = FeatureStore.load(fn);
            delete(fn);
            ex = Experiment();
            ex.FeatureStore = fs;
            ex.Results = struct('Classification', struct(), 'Regression', struct(), 'DimReduction', struct());
        end

    end
end
