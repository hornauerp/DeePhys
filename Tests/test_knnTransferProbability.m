classdef test_knnTransferProbability < matlab.unittest.TestCase
% TEST_KNNTRANSFERPROBABILITY  Inductive transfer of P(inhibitory) to external units.
%
% Run: runtests('test_knnTransferProbability')

    methods (Test)

        function testInverseDistanceWeights(tc)
            % the unclassified reference unit (NaN) is left out
            X = [0; 1; 3; 10];
            p = [0; 1; NaN; 1];
            q = knnTransferProbability(X, p, [0.5; 0.25], 2);
            % 0.5: distances 0.5 and 0.5 -> 0.5; 0.25: weights 4 (p 0) and 4/3 (p 1) -> 0.25
            tc.verifyEqual(q, [0.5; 0.25], 'AbsTol', 1e-12);
        end

        function testNoClassifiedReferenceGivesNaN(tc)
            q = knnTransferProbability([0; 1], [NaN; NaN], [0.5; 2], 3);
            tc.verifyTrue(all(isnan(q)));
        end

        function testShiftedExternalClusterKeepsItsClasses(tc)
            % external units shifted away from the reference data and tight among themselves
            rng(7);
            y = [ones(200, 1); 2 * ones(200, 1)];
            X = randn(400, 3) .* [1 0.3 0.3] + (y - 1.5) .* [4 0 0];
            ye = [ones(80, 1); 2 * ones(20, 1)];
            Xe = randn(100, 3) * 0.3 + (ye - 1.5) .* [4 0 0] + [0 3 0];
            q = knnTransferProbability(X, double(y == 2), Xe, 10);
            tc.verifyGreaterThan(median(q(ye == 2)), 0.9);
            tc.verifyLessThan(median(q(ye == 1)), 0.1);
        end

    end
end
