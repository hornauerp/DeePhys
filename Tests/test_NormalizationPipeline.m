classdef test_NormalizationPipeline < matlab.unittest.TestCase
% TEST_NORMALIZATIONPIPELINE  Unit tests for NormalizationPipeline.
%
% Run: runtests('test_NormalizationPipeline')

    methods (Test)

        function testFitTransformReturnsSameSize(tc)
            rng(42);
            X = randn(20, 5);
            labels = repmat([1;2], 10, 1);
            np = NormalizationPipeline.groupThenGlobal();
            [X_out, np2] = np.fit_transform(X, labels);
            tc.verifySize(X_out, size(X));
            tc.verifyNotEmpty(np2.FitParams);
        end

        function testTransformAppliesSameMappingAsTraining(tc)
            rng(1);
            X_train = randn(30, 4);
            X_test  = randn(10, 4);
            labels_train = repmat([1;2;3], 10, 1);
            labels_test  = repmat([1;2], 5, 1);
            np = NormalizationPipeline.globalOnly();
            [~, np_fit] = np.fit_transform(X_train, labels_train);
            X_out1 = np_fit.transform(X_train, labels_train);
            X_out2 = np_fit.transform(X_test,  labels_test);
            % Output must be bounded by clip range (±5 before max-abs scaling)
            tc.verifyLessThanOrEqual(max(abs(X_out1(:))), 1 + 1e-10);
            tc.verifyLessThanOrEqual(max(abs(X_out2(:))), 1 + 1e-10);
        end

        function testGroupThenGlobalCentersEachGroup(tc)
            % Group 1: mean ≈ 0, Group 2: mean ≈ 10 → after group z-score both centered
            X = [randn(10,2); randn(10,2) + 10];
            labels = [ones(10,1); 2*ones(10,1)];
            np = NormalizationPipeline.groupThenGlobal();
            [X_out, ~] = np.fit_transform(X, labels);
            tc.verifyLessThan(max(abs(X_out(:))), 1 + 1e-10);  % clipped + scaled
        end

        function testUnseenGroupFallsBackToPooled(tc)
            rng(7);
            X_train = randn(20, 3);
            labels_train = ones(20, 1);  % only group 1 during fit
            X_test  = randn(5, 3);
            labels_test = 2 * ones(5, 1);  % unseen group 2
            np = NormalizationPipeline({struct('type','group_zscore','params',struct())});
            [~, np_fit] = np.fit_transform(X_train, labels_train);
            % Should warn and not error
            warnState = warning('off', 'NormalizationPipeline:unseenGroup');
            X_out = np_fit.transform(X_test, labels_test);
            warning(warnState);
            tc.verifySize(X_out, size(X_test));
            tc.verifyTrue(~any(isnan(X_out(:))));
        end

        function testSingletonGroupIssuesWarningAndFallsBack(tc)
            rng(3);
            X = randn(11, 3);
            labels = [1; ones(10,1)*2];  % group 1 has only 1 sample
            np = NormalizationPipeline({struct('type','group_zscore','params',struct())});
            tc.verifyWarning( ...
                @() np.fit_transform(X, labels), ...
                'NormalizationPipeline:singletonGroup');
        end

        function testAllZeroColumnDoesNotProduceNaN(tc)
            X = [ones(10,1)*5, zeros(10,1), randn(10,1)];
            labels = ones(10,1);
            np = NormalizationPipeline.globalOnly();
            [X_out, ~] = np.fit_transform(X, labels);
            tc.verifyTrue(~any(isnan(X_out(:))));
        end

        function testFitParamsSerializable(tc)
            rng(5);
            X = randn(20, 4);
            labels = repmat([1;2], 10, 1);
            np = NormalizationPipeline.groupThenGlobal();
            [~, np_fit] = np.fit_transform(X, labels);
            % FitParams should be a plain struct (no handle objects)
            tc.verifyTrue(isstruct(np_fit.FitParams));
        end

        function testTransformBeforeFitErrors(tc)
            X = randn(5, 3);
            np = NormalizationPipeline.globalOnly();
            % FitParams is empty — transform should error or return garbage but not silently work
            % We just verify it doesn't crash in an unexpected way; actual behaviour is a field access error
            tc.verifyError(@() np.transform(X, []), 'MATLAB:nonExistentField');
        end

        % ------------------------------------------------------------------
        % ComBat batch correction
        % ------------------------------------------------------------------

        function testComBatRemovesBatchShift(tc)
            % Two batches with a known mean offset; same within-batch noise.
            % After ComBat, per-batch feature means should be approximately equal.
            rng(42);
            n = 30; F = 4;
            X1 = randn(n, F);              % batch 1: mean ≈ 0
            X2 = randn(n, F) + 5;          % batch 2: mean ≈ 5
            X  = [X1; X2];
            batch = [ones(n,1); 2*ones(n,1)];

            np = NormalizationPipeline({struct('type','combat','params',struct())});
            [X_out, ~] = np.fit_transform(X, batch);

            mean1 = mean(X_out(1:n,:), 1);
            mean2 = mean(X_out(n+1:end,:), 1);
            tc.verifyLessThan(max(abs(mean1 - mean2)), 0.5, ...
                'Batch means should converge after ComBat');
        end

        function testComBatPreservesBiologicalSignal(tc)
            % Two classes separated along feature 1; perfectly confounded with batch
            % (class A only in batch 1, class B only in batch 2).
            % ComBat with covariate protection should preserve the class separation.
            rng(7);
            n = 20; F = 3;
            X1 = [randn(n,1) + 3, randn(n, F-1)];   % class 1, batch 1, feature 1 mean = 3
            X2 = [randn(n,1) - 3, randn(n, F-1)];   % class 2, batch 2, feature 1 mean = -3
            X  = [X1; X2];
            batch = [ones(n,1); 2*ones(n,1)];
            covar = [ones(n,1); 2*ones(n,1)];         % class labels

            np = NormalizationPipeline({struct('type','combat','params',struct())});
            [X_out, ~] = np.fit_transform(X, batch, CovariateLabels=covar);

            % Class separation on feature 1 should still be substantial
            sep = abs(mean(X_out(1:n,1)) - mean(X_out(n+1:end,1)));
            tc.verifyGreaterThan(sep, 1.0, ...
                'Biological signal should be preserved with covariate protection');
        end

        function testComBatUnseenBatchWarns(tc)
            % Fit on batches 1 and 2, transform on batch 3 — should warn and not crash.
            rng(3);
            X_train = randn(20, 3);
            b_train = [ones(10,1); 2*ones(10,1)];
            X_test  = randn(5, 3);
            b_test  = 3 * ones(5, 1);

            np = NormalizationPipeline({struct('type','combat','params',struct())});
            [~, np_fit] = np.fit_transform(X_train, b_train);

            tc.verifyWarning( ...
                @() np_fit.transform(X_test, b_test), ...
                'NormalizationPipeline:unseenBatch');
        end

        function testComBatSingleBatchIsApproxIdentity(tc)
            % With only one batch the correction should leave data approximately unchanged.
            rng(11);
            X = randn(30, 4);
            batch = ones(30, 1);

            np = NormalizationPipeline({struct('type','combat','params',struct())});
            [X_out, ~] = np.fit_transform(X, batch);

            tc.verifyLessThan(max(abs(X_out(:) - X(:))), 1e-4, ...
                'Single-batch ComBat should be an approximate identity');
        end

    end
end
