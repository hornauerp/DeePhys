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
            [X_fit, np_fit] = np.fit_transform(X_train, labels_train);
            X_out1 = np_fit.transform(X_train, labels_train);
            X_out2 = np_fit.transform(X_test,  labels_test);
            % transform() reproduces the fitted mapping on the training data ...
            tc.verifyEqual(X_out1, X_fit, 'AbsTol', 1e-12);
            tc.verifyLessThanOrEqual(max(abs(X_out1(:))), 1 + 1e-10);
            % ... and applies the same affine map to new data, which may exceed the
            % training range (max-abs scaling is fitted on the training set only).
            tc.verifyTrue(all(isfinite(X_out2(:))));
            ref = np_fit.transform([X_train(1, :); X_test(1, :)], [labels_train(1); labels_test(1)]);
            tc.verifyEqual(ref(2, :), X_out2(1, :), 'AbsTol', 1e-12);
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
            % Class and batch partially confounded (batch 1: 75% class A, batch 2: 25%),
            % class effect +-3 on feature 1, batch offset +4 on all features. With the
            % covariate, ComBat must remove the batch offset but keep the class effect.
            rng(7);
            n = 40; F = 6;
            batch = [ones(n,1); 2*ones(n,1)];
            covar = [ones(30,1); 2*ones(10,1); ones(10,1); 2*ones(30,1)];
            X = randn(2*n, F);
            X(:,1) = X(:,1) + 3 * (covar == 1) - 3 * (covar == 2);
            X(batch == 2, :) = X(batch == 2, :) + 4;

            np = NormalizationPipeline({struct('type','combat','params',struct())});
            [X_out, ~] = np.fit_transform(X, batch, CovariateLabels=covar);

            sep = mean(X_out(covar == 1, 1)) - mean(X_out(covar == 2, 1));
            tc.verifyGreaterThan(sep, 4.5, 'Class effect (6) should be largely preserved');
            within = @(c, b) mean(X_out(covar == c & batch == b, 2:end), 'all');
            tc.verifyLessThan(abs(within(1, 1) - within(1, 2)), 0.5, 'Batch offset should be removed');
        end

        function testComBatConfoundedCovariateWarns(tc)
            % Class identical to batch: the effects cannot be separated.
            rng(7);
            X = [randn(20, 3) + 3; randn(20, 3) - 3];
            batch = [ones(20,1); 2*ones(20,1)];
            np = NormalizationPipeline({struct('type','combat','params',struct())});
            tc.verifyWarning(@() np.fit_transform(X, batch, CovariateLabels=batch), ...
                'NormalizationPipeline:combatConfounded');
        end

        function testComBatTransformIsLabelFree(tc)
            % Labels shape the fit only; transforming the training data again (without
            % labels) must reproduce fit_transform, and labels passed to transform are ignored.
            rng(5);
            X = [randn(30, 4); randn(30, 4) + 2];
            batch = [ones(30,1); 2*ones(30,1)];
            covar = repmat([1; 2], 30, 1);
            np = NormalizationPipeline({struct('type','combat','params',struct())});
            [X_fit, np_fit] = np.fit_transform(X, batch, CovariateLabels=covar);
            tc.verifyEqual(np_fit.transform(X, batch), X_fit, 'AbsTol', 1e-10);
            tc.verifyEqual(np_fit.transform(X, batch, CovariateLabels=3 - covar), X_fit, 'AbsTol', 1e-10);
        end

        function testComBatUnseenBatchIsCorrected(tc)
            % Grouped CV: the held-out batch was not seen during fit. Its offset must still be
            % removed (location/scale estimated from its own samples, without labels).
            rng(3);
            X_train = [randn(20, 5); randn(20, 5) + 3];
            b_train = [ones(20,1); 2*ones(20,1)];
            X_test  = randn(15, 5) + 8;
            b_test  = 3 * ones(15, 1);
            np = NormalizationPipeline({struct('type','combat','params',struct())});
            [X_tr, np_fit] = np.fit_transform(X_train, b_train);
            X_te = np_fit.transform(X_test, b_test);
            tc.verifyLessThan(abs(mean(X_te, 'all') - mean(X_tr, 'all')), 0.5);
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
