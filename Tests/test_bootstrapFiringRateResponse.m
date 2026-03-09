classdef test_bootstrapFiringRateResponse < matlab.unittest.TestCase
% TEST_BOOTSTRAPFIRINGRATERESPONSE  Unit tests for Functions/bootstrapFiringRateResponse.m
%
% Run: runtests('test_bootstrapFiringRateResponse')

    methods (Test)

        function testOutputStructHasRequiredFields(tc)
            pre  = rand(4, 20);
            post = rand(4, 20);
            r    = bootstrapFiringRateResponse(pre, post, 100, 0.05);
            tc.verifyTrue(isfield(r, 'increase'));
            tc.verifyTrue(isfield(r, 'decrease'));
            tc.verifyTrue(isfield(r, 'unchanged'));
        end

        function testAllUnitsCategorizingExactlyOnce(tc)
            rng(42);
            n_units = 10; n_bins = 50;
            pre  = rand(n_units, n_bins);
            post = rand(n_units, n_bins);
            r    = bootstrapFiringRateResponse(pre, post, 500, 0.05);
            all_idx = sort([r.increase, r.decrease, r.unchanged]);
            tc.verifyEqual(all_idx, 1:n_units, ...
                'Every unit must appear in exactly one category');
        end

        function testClearIncreaseDetected(tc)
            % post = 10x pre — huge step-up, all units must be flagged as increased
            rng(0);
            n_units = 5; n_bins = 100;
            pre  = zeros(n_units, n_bins);
            post = 10 * ones(n_units, n_bins);
            r    = bootstrapFiringRateResponse(pre, post, 1000, 0.01);
            tc.verifyEqual(sort(r.increase), 1:n_units);
            tc.verifyEmpty(r.decrease);
        end

        function testClearDecreaseDetected(tc)
            % post = 0, pre = 10 — huge drop, all units must be flagged as decreased
            rng(1);
            n_units = 5; n_bins = 100;
            pre  = 10 * ones(n_units, n_bins);
            post = zeros(n_units, n_bins);
            r    = bootstrapFiringRateResponse(pre, post, 1000, 0.01);
            tc.verifyEqual(sort(r.decrease), 1:n_units);
            tc.verifyEmpty(r.increase);
        end

        function testIdenticalMatricesAllUnchanged(tc)
            % Empirical diff is exactly 0 — must sit inside the null interval
            rng(7);
            n_units = 8; n_bins = 40;
            mat = rand(n_units, n_bins);
            r   = bootstrapFiringRateResponse(mat, mat, 500, 0.01);
            tc.verifyEmpty(r.increase);
            tc.verifyEmpty(r.decrease);
            tc.verifyEqual(sort(r.unchanged), 1:n_units);
        end

        function testDimensionMismatchThrowsError(tc)
            pre  = rand(3, 10);
            post = rand(3, 15);   % different bin count
            tc.verifyError(@() bootstrapFiringRateResponse(pre, post), ...
                'MATLAB:assert:failed');
        end

        function testSingleUnitWorks(tc)
            rng(3);
            pre  = rand(1, 30);
            post = rand(1, 30);
            r    = bootstrapFiringRateResponse(pre, post, 200, 0.05);
            all_idx = sort([r.increase, r.decrease, r.unchanged]);
            tc.verifyEqual(all_idx, 1);
        end

        function testDefaultArgumentsAccepted(tc)
            % Calling with only two arguments should not error (uses n_iter=1000, alpha=0.01)
            pre  = rand(3, 20);
            post = rand(3, 20);
            tc.verifyWarningFree(@() bootstrapFiringRateResponse(pre, post));
        end

        function testIndexReturnedNotBool(tc)
            % increase/decrease/unchanged must contain integer indices, not logicals
            rng(9);
            pre  = rand(6, 30);
            post = rand(6, 30);
            r    = bootstrapFiringRateResponse(pre, post, 300, 0.05);
            tc.verifyTrue(isnumeric(r.increase));
            tc.verifyTrue(isnumeric(r.decrease));
            tc.verifyTrue(isnumeric(r.unchanged));
        end

    end
end
