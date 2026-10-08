classdef test_spikeRhythmicity < matlab.unittest.TestCase
% TEST_SPIKERHYTHMICITY  Rhythm detection on simulated spike trains with known truth.
%
% Run: runtests('test_spikeRhythmicity')

    methods (Test)

        function testSineModulatedRhythmDetected(tc)
            rng(11);
            st = inhomPoisson(@(t) 5 * (1 + 0.5 * sin(2 * pi * 8 * t)), 600, 10);
            r = spikeRhythmicity(st, 600);
            tc.verifyTrue(r.Rhythmic);
            tc.verifyEqual(r.Freq, 8, 'AbsTol', 0.5);
        end

        function testPoissonNotRhythmic(tc)
            rng(12);
            calls = false(1, 20);
            for i = 1:20
                st = sort(rand(3000, 1) * 600);
                calls(i) = spikeRhythmicity(st, 600).Rhythmic;
            end
            tc.verifyLessThanOrEqual(sum(calls), 2);
            r = spikeRhythmicity(sort(rand(3000, 1) * 600), 600);
            tc.verifyEqual(r.LV, 1, 'AbsTol', 0.1);                % Poisson: LV ~ 1
        end

        function testRegularFiringHasLowLV(tc)
            rng(13);
            st = (0.1:0.1:600)' + 0.005 * randn(6000, 1);
            r = spikeRhythmicity(st, 600);
            tc.verifyLessThan(r.LV, 0.2);
            tc.verifyTrue(r.Rhythmic);
            tc.verifyEqual(r.Freq, 10, 'AbsTol', 0.5);
        end

        function testTooFewSpikes(tc)
            r = spikeRhythmicity([1; 2; 3], 600);
            tc.verifyFalse(r.Rhythmic);
            tc.verifyTrue(isnan(r.Freq));
        end

    end
end

function st = inhomPoisson(rate_fn, T, max_rate)
% thinning: candidates at max_rate, kept with probability rate / max_rate
    cand = sort(rand(poissrnd(max_rate * T), 1) * T);
    st = cand(rand(numel(cand), 1) < rate_fn(cand) / max_rate);
end
