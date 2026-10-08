classdef test_burstEI < matlab.unittest.TestCase
% TEST_BURSTEI  burstEIMetrics, burstUnitTiming and burstRhythm on synthetic networks with known
% E / I behaviour around the bursts.
%
% Run: runtests('test_burstEI')

    methods (Test)

        function testIEDropAndILeadBeyondRateMatchedNull(tc)
            [st, su, lab, b, dur] = test_burstEI.network(1);
            R = burstEIMetrics(st, su, lab, b, dur, NNull=100);
            tc.verifyTrue(R.Valid);
            m = test_burstEI.metric(R);
            tc.verifyLessThan(m.IEContrastIn.Value, m.IEContrastIn.NullLow);      % I/E lower in bursts
            tc.verifyLessThan(m.IEContrastIn.PLower, 0.05);
            tc.verifyGreaterThan(m.IFoldPre.Value, 2);                            % I fires before onset
            tc.verifyLessThan(m.IFoldIn.Value, 0.5);                              % and is silent inside
            tc.verifyLessThan(m.ILeadTime.Value, 0);                             % I fires uniformly in
            tc.verifyGreaterThanOrEqual(m.ILeadTime.Value, -0.2);                % the 0.2 s before onset
            tc.verifyGreaterThan(m.IEMaxTime.Value, -0.25);
            tc.verifyLessThan(m.IEMaxTime.Value, 0);
            tc.verifyTrue(isnan(m.AntiPhase.Value));                             % no superbursts
            % traces: I/E of the onset-aligned average dips at the burst
            t = R.Traces.Onset; k = t.Time > 0.02 & t.Time < 0.12; p = t.Time > -0.2 & t.Time < -0.02;
            tc.verifyLessThan(mean(t.IE(k)), 0.1 * mean(t.IE(p)));
        end

        function testIAfterPeak(tc)
            % I rebounds 50-100 ms after each burst offset (0.125-0.175 s after the peak); E not
            [st, su, lab, b, dur] = test_burstEI.network(10, Rebound=true);
            m = test_burstEI.metric(burstEIMetrics(st, su, lab, b, dur, NNull=100));
            tc.verifyGreaterThan(m.IAfterPeak.Value, m.IAfterPeak.NullHigh);
            tc.verifyLessThan(abs(m.IAfterPeakTime.Value - 0.15), 0.035);
            tc.verifyLessThan(m.EAfterPeak.Value, m.IAfterPeak.Value / 3);
            % without the rebound I only decays after the burst
            [st, su, lab, b, dur] = test_burstEI.network(10);
            m0 = test_burstEI.metric(burstEIMetrics(st, su, lab, b, dur, NNull=20));
            tc.verifyLessThan(m0.IAfterPeak.Value, m.IAfterPeak.Value / 3);
        end

        function testUninformativeLabelsStayInsideNull(tc)
            % "I" = random E units: the label carries nothing beyond firing rate
            [st, su, lab, b, dur] = test_burstEI.network(2);
            rs = RandStream('twister', 'Seed', 5);
            e = find(lab == 1); fake = nan(size(lab)); fake(e) = 1;
            fake(e(randperm(rs, numel(e), 20))) = 2;
            R = burstEIMetrics(st, su, fake, b, dur, NNull=100);
            m = test_burstEI.metric(R);
            tc.verifyGreaterThan(m.IEContrastIn.PLower, 0.025);
            tc.verifyGreaterThan(m.IEContrastIn.PUpper, 0.025);
        end

        function testTooFewUnitsIsInvalid(tc)
            [st, su, lab, b, dur] = test_burstEI.network(3);
            lab(lab == 2) = NaN; lab(find(isnan(lab), 2)) = 2;
            R = burstEIMetrics(st, su, lab, b, dur, NNull=10);
            tc.verifyFalse(R.Valid);
            tc.verifyTrue(all(isnan(R.Metrics.Value)));
        end

        function testUnitTiming(tc)
            [st, su, lab, b, dur] = test_burstEI.network(4);
            T = burstUnitTiming(st, su, b, dur);
            tc.verifyEqual(height(T), numel(lab));
            tc.verifyLessThan(median(T.Latency(lab == 2), 'omitnan'), -0.05);
            tc.verifyGreaterThan(median(T.Latency(lab == 1), 'omitnan'), -0.02);
            tc.verifyGreaterThan(median(T.Participation(lab == 1)), 0.8);
            tc.verifyLessThan(median(T.Participation(lab == 2)), 0.3);
            tc.verifyGreaterThan(median(T.PreRatio(lab == 2)), 3);
            tc.verifyLessThan(median(T.PreRatio(lab == 1)), 1.5);
        end

        function testAntiPhaseInsideSuperbursts(tc)
            [st, su, lab, b, dur] = test_burstEI.network(5, Superbursts=true);
            R = burstEIMetrics(st, su, lab, b, dur, NNull=100);
            m = test_burstEI.metric(R);
            tc.verifyGreaterThanOrEqual(R.NInnerBursts, 5);
            tc.verifyGreaterThan(m.AntiPhase.Value, m.AntiPhase.NullHigh);
            tc.verifyLessThan(m.AntiPhase.PUpper, 0.05);
        end

        function testLockedSubPeaks(tc)
            % 10 Hz sub-peaks at fixed latencies from the burst peak; a smooth envelope does not
            % produce reliable fast structure
            [st, su, b] = test_burstEI.rhythmic(6, 10, "locked");
            R = burstRhythm(st, su, struct('All', true(1, 60)), b);
            L = R.Locked;
            tc.verifyLessThan(abs(L.Freq - 10), 1.5);
            tc.verifyGreaterThan(L.Reliability, 0.8);
            tc.verifyLessThan(L.P, 0.05);
            [st, su, b] = test_burstEI.rhythmic(7, 10, "none");
            L = burstRhythm(st, su, struct('All', true(1, 60)), b).Locked;
            tc.verifyTrue(L.P > 0.05 || L.Reliability < 0.5);
        end

        function testSharedRhythmNotLocked(tc)
            % a 15 Hz rhythm with a random phase in every burst: shared between disjoint units,
            % invisible in the burst average
            [st, su, b] = test_burstEI.rhythmic(8, 15, "random");
            R = burstRhythm(st, su, struct('All', true(1, 60)), b, NShift=20);
            S = R.Shared;
            tc.verifyLessThan(abs(S.PeakFreq - 15), 3);
            tc.verifyLessThan(S.P, 0.1);
            tc.verifyGreaterThan(S.Coh15_30, 0.1);
            [st, su, b] = test_burstEI.rhythmic(9, 15, "none");
            S0 = burstRhythm(st, su, struct('All', true(1, 60)), b, NShift=20).Shared;
            tc.verifyLessThan(abs(S0.Coh15_30), 0.5 * S.Coh15_30);
        end
    end

    methods (Static)
        function m = metric(R)
            m = struct();
            for k = 1:height(R.Metrics), m.(R.Metrics.Metric(k)) = table2struct(R.Metrics(k, 2:end)); end
        end

        function [st, su, lab, b, dur] = network(seed, opts)
            % 60 E units (bursts at 40 Hz), 20 I units (20 Hz in the 0.2 s before each onset,
            % nearly silent in bursts); background rates 0.5-4 Hz in both, so rates overlap.
            % Superbursts: bursts in groups of 5, 0.5 s apart, I firing in the gaps instead.
            % Rebound: I also fires at 15 Hz 50-100 ms after each burst offset.
            arguments
                seed (1,1) double
                opts.Superbursts (1,1) logical = false
                opts.Rebound (1,1) logical = false
            end
            rs = RandStream('twister', 'Seed', seed);
            dur = 300; dt = 1e-3; ng = dur / dt; nE = 60; nI = 20;
            if opts.Superbursts
                c = 5:6:290; on = reshape(c + (0:4)' * 0.5, 1, []); sbi = repelem(1:numel(c), 5);
            else
                on = cumsum(3 + rs.rand(1, 80)); on = on(on < dur - 2); sbi = zeros(size(on));
            end
            off = on + 0.15;
            inB = false(ng, 1); pre = false(ng, 1); gap = false(ng, 1); reb = false(ng, 1);
            for k = 1:numel(on)
                inB(round(on(k) / dt) + 1:round(off(k) / dt)) = true;
                if opts.Rebound, reb(round((off(k) + 0.05) / dt) + 1:round((off(k) + 0.1) / dt)) = true; end
                if ~opts.Superbursts || sbi(k) ~= sbi(max(k - 1, 1)) || k == 1
                    pre(max(1, round((on(k) - 0.2) / dt) + 1):round(on(k) / dt)) = true;
                else
                    gap(round(off(k - 1) / dt) + 1:round(on(k) / dt)) = true;
                end
            end
            st = []; su = [];
            for u = 1:nE + nI
                bg = 0.5 + 3.5 * rs.rand;
                r = bg * ones(ng, 1);
                if u <= nE, r(inB) = 40;
                else, r(pre | gap) = 20; r(reb) = 15; r(inB) = 0.2; end
                t = (find(rs.rand(ng, 1) < r * dt) - rs.rand) * dt;
                st = [st; t]; su = [su; u * ones(numel(t), 1)]; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
            lab = [ones(1, nE), 2 * ones(1, nI)];
            b = struct('Onset', on, 'Offset', off, 'T_peak', on + 0.075, 'Superburst', sbi);
        end

        function [st, su, b] = rhythmic(seed, f, mode)
            % 60 units, 100 bursts with a Gaussian envelope (SD 0.1 s, peak 40 Hz per unit);
            % mode "locked": x (1 + 0.9 cos(2 pi f (t - peak))); "random": a random phase per
            % burst; "none": envelope only
            rs = RandStream('twister', 'Seed', seed);
            nu = 60; nb = 100; pk = 2 + (0:nb - 1) * 3 + 0.5 * rs.rand(1, nb);
            dur = pk(end) + 2; dt = 5e-4; t = (0:dt:dur)';
            env = zeros(size(t)); ph = zeros(size(t));
            for k = 1:nb
                i = abs(t - pk(k)) < 0.5;
                env(i) = env(i) + exp(-0.5 * ((t(i) - pk(k)) / 0.1).^2);
                switch mode
                    case "locked", ph(i) = 2 * pi * f * (t(i) - pk(k));
                    case "random", ph(i) = 2 * pi * f * (t(i) - pk(k)) + 2 * pi * rs.rand;
                end
            end
            r = 0.5 + 40 * env;
            if mode ~= "none", r = 0.5 + 40 * env .* (1 + 0.9 * cos(ph)); end
            st = []; su = [];
            for u = 1:nu
                s = t(rs.rand(numel(t), 1) < r * dt) + dt * rs.rand;
                st = [st; s]; su = [su; u * ones(numel(s), 1)]; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
            b = struct('Onset', pk - 0.2, 'Offset', pk + 0.2, 'T_peak', pk);
        end
    end
end
