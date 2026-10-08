classdef test_acgSubtypeLegacy < matlab.unittest.TestCase
% TEST_ACGSUBTYPELEGACY  Behaviour of the legacy ACG subtype rules on clear cases.
%
% Run: runtests('test_acgSubtypeLegacy')

    methods (Test)

        function testRegularFiringIsOscillatory(tc)
            % spikes every 125 ms (8 Hz) with 3 ms jitter for 600 s
            rng(1);
            st = (0.5:0.125:600)' + 0.003 * randn(numel(0.5:0.125:600), 1);
            [cls, info] = acgSubtypeLegacy(acgOf(st), 0.002);
            tc.verifyEqual(cls, "oscillatory");
            tc.verifyEqual(info.PeakFreqACG, 8, 'AbsTol', 0.3);
        end

        function testEmptyACGIsNoisy(tc)
            tc.verifyEqual(acgSubtypeLegacy(zeros(1, 2001), 0.002), "noisy");
        end

        function testOneRowPerUnit(tc)
            rng(2);
            A = [acgOf(sort(rand(3000, 1) * 600)); acgOf(sort(rand(3000, 1) * 600))];
            [cls, info] = acgSubtypeLegacy(A, 0.002);
            tc.verifySize(cls, [2 1]);
            tc.verifySize(info, [2 1]);
        end

    end
end

function a = acgOf(st)
    a = CellTypeClassifier.computeACG(st, 0.002, 2, 10000)';
    a = a / max(max(a), 1);
end
