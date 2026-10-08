classdef test_burstShape < matlab.unittest.TestCase
% TEST_BURSTSHAPE  Burst onset / offset / shape at a fixed fraction of each burst's own amplitude.
%
% Run: runtests('test_burstShape')

    methods (Test)

        function testLinearRampsGiveExpectedTimes(tc)
            % triangular bursts: linear rise over 100 ms, fall over 200 ms; at Level L the onset
            % is (1 - L) x 100 ms before the peak, the offset (1 - L) x 200 ms after it
            [st, pk] = test_burstShape.triangularBursts(1);
            for L = [0.25 0.5]
                S = burstShape(st, pk - 0.15, pk + 0.25, 60, Level=L);
                tc.verifyEqual(median(S.Peak - pk(:)), 0, 'AbsTol', 0.01);   % smoothing shifts the peak of an asymmetric ramp
                tc.verifyEqual(median(S.Rise), (1 - L) * 0.1, 'AbsTol', 0.01);
                tc.verifyEqual(median(S.Fall), (1 - L) * 0.2, 'AbsTol', 0.015);
            end
            S = burstShape(st, pk - 0.15, pk + 0.25, 60);
            tc.verifyEqual(median(S.FWHM), 0.15, 'AbsTol', 0.015);
        end

        function testIndependentOfDetectorEdges(tc)
            % the same bursts given with wide, narrow or shifted detector spans: identical shape
            [st, pk] = test_burstShape.triangularBursts(2);
            A = burstShape(st, pk - 0.3, pk + 0.4, 60);
            B = burstShape(st, pk - 0.02, pk + 0.03, 60);
            C = burstShape(st, pk - 0.1, pk + 0.1, 60);
            tc.verifyEqual(B.Onset, A.Onset, 'AbsTol', 0.003);
            tc.verifyEqual(C.Offset, A.Offset, 'AbsTol', 0.003);
            tc.verifyEqual(B.Duration, A.Duration, 'AbsTol', 0.005);
        end

        function testSubEventsDoNotCutTheBurst(tc)
            % two triangular sub-events 150 ms apart (amplitudes 1500 and 2000) inside one span:
            % onset on the first one's ramp, offset on the second one's decay; FWHM stays local
            rs = RandStream('twister', 'Seed', 4);
            T = 30; dt = 0.0005; t = (0:dt:T - dt)'; rate = 20 * ones(size(t)); pk = 2:3:T - 2;
            for p = pk
                for q = [p, p + 0.15]
                    a = 1500 + 500 * (q > p);
                    up = t >= q - 0.05 & t < q; rate(up) = rate(up) + a * (t(up) - (q - 0.05)) / 0.05;
                    dn = t >= q & t < q + 0.05; rate(dn) = rate(dn) + a * (1 - (t(dn) - q) / 0.05);
                end
            end
            hit = rand(rs, size(t)) < rate * dt; st = sort(t(hit) + dt * rand(rs, sum(hit), 1));
            S = burstShape(st, pk - 0.06, pk + 0.21, T);
            % 25 % of 2000 is reached at a third of the first (1500) ramp; 75 % down the last decay
            tc.verifyEqual(median(S.Onset - (pk(:) - 0.05)), 0.05 * 500 / 1500, 'AbsTol', 0.01);
            tc.verifyEqual(median(S.Offset - (pk(:) + 0.2)), -0.05 * 0.25, 'AbsTol', 0.01);
            tc.verifyLessThan(median(S.FWHM), 0.08);
        end

        function testNestedBurstNotCutAtFirstDip(tc)
            % five sub-bursts 100 ms apart with silence between them and no background (the first
            % rate minimum after the peak is a dip inside the burst): onset / offset still span
            % all five
            rs = RandStream('twister', 'Seed', 5);
            T = 40; dt = 0.0005; t = (0:dt:T - dt)'; rate = zeros(size(t)); pk = 3:5:T - 3;
            for p = pk
                for q = p + (0:4) * 0.1
                    in = t >= q & t < q + 0.02; rate(in) = rate(in) + 2000 * (1 + 0.5 * (q == p));
                end
            end
            hit = rand(rs, size(t)) < rate * dt; st = sort(t(hit) + dt * rand(rs, sum(hit), 1));
            S = burstShape(st, pk - 0.01, pk + 0.43, T);
            tc.verifyEqual(median(S.Onset - pk(:)), 0, 'AbsTol', 0.02);
            tc.verifyEqual(median(S.Offset - (pk(:) + 0.42)), 0, 'AbsTol', 0.02);
        end

        function testOrderAndEmptyInput(tc)
            [st, pk] = test_burstShape.triangularBursts(3);
            o = numel(pk):-1:1;                               % reversed input order
            A = burstShape(st, pk - 0.15, pk + 0.25, 60);
            B = burstShape(st, pk(o) - 0.15, pk(o) + 0.25, 60);
            tc.verifyEqual(B.Peak, A.Peak(o), 'AbsTol', 1e-12);
            E = burstShape(st, zeros(1, 0), zeros(1, 0), 60);
            tc.verifyEqual(height(E), 0);
        end

    end

    methods (Static)
        function [st, pk] = triangularBursts(seed)
            % Poisson background (20 spikes/s) plus bursts with a triangular rate profile
            % (0 -> 2000 spikes/s over 100 ms, back to 0 over 200 ms), every 3 s
            rs = RandStream('twister', 'Seed', seed);
            T = 60; dt = 0.0005; t = (0:dt:T - dt)';
            pk = 1.5:3:T - 1;
            rate = 20 * ones(size(t));
            for p = pk
                up = t >= p - 0.1 & t < p;  rate(up) = rate(up) + 2000 * (t(up) - (p - 0.1)) / 0.1;
                dn = t >= p & t < p + 0.2;  rate(dn) = rate(dn) + 2000 * (1 - (t(dn) - p) / 0.2);
            end
            hit = rand(rs, size(t)) < rate * dt;
            st = t(hit) + dt * rand(rs, sum(hit), 1);
            st = sort(st);
        end
    end
end
