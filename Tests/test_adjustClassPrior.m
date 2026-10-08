classdef test_adjustClassPrior < matlab.unittest.TestCase
% TEST_ADJUSTCLASSPRIOR  EM prior re-estimation recovers class fractions under prior shift.
%
% Run: runtests('test_adjustClassPrior')

    methods (Test)

        function testRecoversShiftedPrior(tc)
            % two Gaussian classes (d' = 2), classifier trained at 50 %, applied at 10 % and 80 %
            rng(21);
            post = @(x) 1 ./ (1 + exp(-2 * x));                 % exact posterior at prior 0.5
            for q_true = [0.1 0.8]
                n = 4000; y = rand(n, 1) < q_true;
                x = randn(n, 1) + (2 * y - 1);
                q = adjustClassPrior(post(x), 0.5);
                tc.verifyEqual(q, q_true, 'AbsTol', 0.03);
                tc.verifyGreaterThan(abs(mean(post(x) > 0.5) - q_true), abs(q - q_true));   % beats counting
            end
        end

        function testUnchangedPriorIsFixedPoint(tc)
            rng(22);
            p = rand(500, 1);
            [q, post] = adjustClassPrior(p, mean(p), Init=mean(p));
            tc.verifyEqual(q, mean(p), 'AbsTol', 1e-6);
            tc.verifyEqual(post, p, 'AbsTol', 1e-6);
        end

        function testEmptyInput(tc)
            tc.verifyTrue(isnan(adjustClassPrior(zeros(0, 1), 0.5)));
        end

    end
end
