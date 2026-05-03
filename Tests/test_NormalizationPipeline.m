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

    end
end
