classdef test_FeatureFixes < matlab.unittest.TestCase
% TEST_FEATUREFIXES  Feature definitions checked against synthetic data with a known answer.
%
% Run: runtests('test_FeatureFixes')

    methods (Test)

        function testSTTCSurrogatesDetectBurstComodulation(tc)
            % six units sharing 150 bursts (spikes jittered by 25 ms) plus 1 Hz background:
            % with surrogates jittered by less than the coincidence window nothing was significant
            rng(1); Ts = 300; n = 6; bt = sort(rand(150, 1) * Ts);
            [st, su, ust] = test_FeatureFixes.trains(n, @() [repelem(bt, 4) + 0.025 * randn(600, 1); rand(Ts, 1) * Ts], Ts);
            sp = struct('MaxLag', 0.05, 'SurrogateMethod', 'jitter', 'N_Surrogates', 50, 'Percentile', 95);
            c = computeConnectivitySTTC(ust, st, su, n, Ts, 20000, sp);
            tc.verifyGreaterThanOrEqual(nnz(triu(c.bu, 1)), 13);   % of 15 pairs
        end

        function testSTTCIndependentUnitsStaySparse(tc)
            rng(2); Ts = 300; n = 12;
            [st, su, ust] = test_FeatureFixes.trains(n, @() rand(poissrnd(5 * Ts), 1) * Ts, Ts);
            sp = struct('MaxLag', 0.05, 'SurrogateMethod', 'jitter', 'N_Surrogates', 50, 'Percentile', 95);
            c = computeConnectivitySTTC(ust, st, su, n, Ts, 20000, sp);
            tc.verifyLessThan(nnz(triu(c.bu, 1)) / (n * (n - 1) / 2), 0.2);
        end

        function testHalfWidthIgnoresLaterDip(tc)
            fs = 20000; t = (0:81)' / fs * 1000;
            wf  = -exp(-((t - 1.5) / 0.12).^2) + 0.35 * exp(-((t - 2.0) / 0.25).^2);
            wf2 = wf - 0.6 * exp(-((t - 3.7) / 0.1).^2);
            f = computeWaveformFeatures([wf wf2] ./ max(abs([wf wf2])), fs);
            tc.verifyEqual(f{2}.HalfWidth, f{1}.HalfWidth, 'AbsTol', 0.02);
            tc.verifyLessThan(f{1}.HalfWidth, 0.4);
        end

        function testTroughToPeakOfBroadUnit(tc)
            fs = 20000; t = (0:81)' / fs * 1000;
            wf = -exp(-((t - 1.5) / 0.15).^2) + 0.3 * exp(-((t - 2.8) / 0.35).^2);   % peak 1.3 ms after the trough
            f = computeWaveformFeatures(wf / max(abs(wf)), fs);
            tc.verifyEqual(f{1}.T2Pdelay, 1.3, 'AbsTol', 0.1);
        end

        function testWaveformAUCIndependentOfSamplingRate(tc)
            g = @(t) -exp(-((t - 1.5) / 0.12).^2) + 0.35 * exp(-((t - 2.0) / 0.25).^2);
            w20 = g((0:81)' / 20); w30 = g((0:122)' / 30);
            a = computeWaveformFeatures(w20 / max(abs(w20)), 20000);
            b = computeWaveformFeatures(w30 / max(abs(w30)), 30000);
            tc.verifyEqual(b{1}.AUC_trough, a{1}.AUC_trough, 'RelTol', 0.05);
        end

        function testSpatialMetricsUnderRandomLabels(tc)
            rng(3);
            [ex, ey] = meshgrid(0:17.5:17.5*219, 0:17.5:17.5*119); ec = [ex(:) ey(:)];
            I = nan(15, 1); m20 = I; m40 = I;
            for r = 1:15
                units = struct('ReferenceElectrode', num2cell(randperm(size(ec, 1), 150)));
                lab = 1 + (rand(1, 150) < 0.2);
                tb = computeSpatialEIBalance(ec, units, lab, struct());
                I(r) = tb.SpatialEIMoransI;
                nw = computeSpatialFeatures(ec, units, lab, [], struct());   m20(r) = nw.SpatialMixingIndex;
                nw = computeSpatialFeatures(ec, units, 1 + (rand(1, 150) < 0.4), [], struct()); m40(r) = nw.SpatialMixingIndex;
            end
            tc.verifyLessThan(abs(mean(I)), 0.03);            % the smoothed grid gave ~0.11
            tc.verifyEqual(mean(m20), 1, 'AbsTol', 0.1);      % composition no longer sets the index
            tc.verifyEqual(mean(m40), 1, 'AbsTol', 0.1);
        end

        function testRichClubColumnsFixed(tc)
            C = zeros(6); C(1, 2) = 1; C(2, 3) = 1; C(3, 1) = 1; C(4, 5) = 1; C(5, 6) = 1; C(6, 4) = 1;
            conn.CCG.bd = C;
            nw = computeGraphFeatures(conn, struct('SmallWorldnessIterations', 2));
            names = string(nw.Properties.VariableNames);
            tc.verifyEqual(sum(startsWith(names, "RichClub")), 5);
        end

        function testRegularityMagnitudeIndependentOfDuration(tc)
            rng(4); m = zeros(1, 2); T = [300 1200];
            for k = 1:2
                x = histcounts(rand(5 * T(k), 1) * T(k), 0:0.1:T(k));
                [~, m(k)] = computeRegularity(x / max(x), 0.1, 20);
            end
            tc.verifyEqual(m(2) / m(1), 1, 'AbsTol', 0.35);   % the raw FFT peak grew ~2x
        end

        function testBootstrapEmpiricalPValueNeverZero(tc)
            % sparse counts give non-Normal nulls (empirical branch at alpha < 1e-3); the strong
            % responder exceeds every permutation, which gave p = 0 without the +1 correction
            rng(5);
            pre = poissrnd(0.05, 40, 30); post = poissrnd(0.05, 40, 30); post(1, :) = poissrnd(5, 1, 30);
            w = warning('off', 'bootstrapFiringRateResponse:nonNormalNull'); c = onCleanup(@() warning(w));
            r = bootstrapFiringRateResponse(pre, post, 200, 1e-4, 'NetworkCorrection', false);
            tc.verifyGreaterThan(min(r.p_values), 0);
            tc.verifyEqual(r.increase(1), 1);
        end

        function testBurstIETraceWithoutEActivity(tc)
            % E fires only inside bursts, I also before them: the pre-onset I/E is undefined
            rng(6); nE = 10; nI = 4; T = 200; on = (5:5:190)';
            st = []; su = [];
            for u = 1:nE + nI
                s = repelem(on, 3) + 0.05 + 0.02 * randn(numel(on) * 3, 1);
                if u > nE, s = [s; on - 0.25 + 0.01 * randn(size(on))]; end %#ok<AGROW>
                st = [st; s]; su = [su; u * ones(size(s))]; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
            B = struct('Onset', on', 'T_peak', on' + 0.05, 'Offset', on' + 0.15, 'T_start', on', 'T_end', on' + 0.15);
            R = burstEIMetrics(st, su, [ones(1, nE) 2 * ones(1, nI)], B, T, NNull=10);
            ie = R.Traces.Onset.IE;
            tc.verifyLessThan(max(ie, [], 'omitnan'), 1e6);   % I ./ eps gave ~1e16
        end

    end

    methods (Static)

        function [st, su, ust] = trains(n, gen, Ts)
            st = []; su = [];
            for u = 1:n
                s = gen(); s = sort(s(s > 0 & s < Ts));
                st = [st; s]; su = [su; u * ones(size(s))]; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
            ust = arrayfun(@(u) st(su == u), 1:n, 'un', 0);
        end

    end
end
