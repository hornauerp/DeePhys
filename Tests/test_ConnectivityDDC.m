classdef test_ConnectivityDDC < matlab.unittest.TestCase
% TEST_CONNECTIVITYDDC  Sparse DDC equals the former dense computation; weighted and
% binary graphs and their graph features.
%
% Run: runtests('test_ConnectivityDDC')

    methods (Test)

        function testSparseMatchesDense(tc)
            [st, su, n, dur] = test_ConnectivityDDC.makeSpikes(13);
            for thr = [0, 0.5, -0.2]
                p = struct('BinSize', 0.001, 'Threshold', thr, 'N_Surrogates', 12, 'Percentile', 95);
                rng(21); c = computeConnectivityDDC(st, su, n, dur, p);
                rng(21); [W, T, B] = test_ConnectivityDDC.denseReference(st, su, n, dur, p);
                off = ~eye(n);
                tc.verifyEqual(c.wd(off), W(off), 'RelTol', 1e-9, 'AbsTol', 1e-12);
                tc.verifyEqual(c.threshold, T, 'RelTol', 1e-9, 'AbsTol', 1e-12);
                tc.verifyEqual(c.bd, B);
                tc.verifyEqual(c.wd_sig, c.wd .* c.bd);
            end
        end

        function testDrivenLinkIsDetected(tc)
            % unit 2 fires 1 ms (one bin) after most spikes of unit 1; DDC relates the activity
            % derivative to activity at zero lag, so it sees lags of about one bin
            rng(3);
            dur = 120; t1 = sort(rand(1500, 1) * dur);
            t2 = sort([t1(rand(1500, 1) < 0.7) + 0.001; rand(300, 1) * dur]);
            t3 = sort(rand(1200, 1) * dur);
            st = [t1; t2; t3]; su = [ones(numel(t1), 1); 2 * ones(numel(t2), 1); 3 * ones(numel(t3), 1)];
            [st, o] = sort(st); su = su(o);
            c = computeConnectivityDDC(st, su, 3, dur, struct('BinSize', 0.001, 'Threshold', 0, ...
                'N_Surrogates', 30, 'Percentile', 95));
            tc.verifyTrue(c.bd(1, 2) == 1 || c.bd(2, 1) == 1);
            tc.verifyNotEqual(c.wd_sig(1, 2) + c.wd_sig(2, 1), 0);
        end

        function testSparseGraphDoesNotHang(tc)
            % A graph whose only links are one reciprocal pair: the BCT rewiring randomiser
            % (small-worldness null model) looped forever on it.
            C = zeros(5); C(3, 4) = 1; C(4, 3) = 1;
            nw = computeGraphFeatures(struct('X', struct('bd', C)), struct('SmallWorldnessIterations', 3));
            tc.verifyTrue(isnan(nw.SmallWorldness_X));
        end

        function testWeightedAndBinaryGraphFeatures(tc)
            [st, su, n, dur] = test_ConnectivityDDC.makeSpikes(5);
            c = computeConnectivityDDC(st, su, n, dur, struct('BinSize', 0.001, 'Threshold', 0, ...
                'N_Surrogates', 20, 'Percentile', 80));
            conn = struct('DDC', c);
            [nw_b, unit_b] = computeGraphFeatures(conn, struct('SmallWorldnessIterations', 2));
            tc.verifyFalse(any(endsWith(string(nw_b.Properties.VariableNames), "_w")));
            [nw, unit] = computeGraphFeatures(conn, struct('SmallWorldnessIterations', 2, ...
                'Types', ["binary", "weighted"]));
            names = string(nw.Properties.VariableNames);
            tc.verifyTrue(all(ismember(["Density_DDC", "MeanStrength_DDC_w", "GlobalEfficiency_DDC_w", ...
                "Modularity_DDC_w"], names)));
            tc.verifyTrue(all(ismember(["Strength_DDC_w", "ClusteringCoefficient_DDC_w", ...
                "Betweenness_DDC_w"], string(unit.Properties.VariableNames))));
            tc.verifyEqual(height(unit), n);
            tc.verifyTrue(all(isfinite(unit.Strength_DDC_w)));
            [nw_w, ~] = computeGraphFeatures(conn, struct('Types', "weighted"));
            tc.verifyTrue(all(endsWith(string(nw_w.Properties.VariableNames), "_w")));
            tc.verifyEqual(width(unit), width(unit_b) + 5);
        end

    end

    methods (Static)

        function [st, su, n, dur] = makeSpikes(n)
            % n units, Poisson 4-12 Hz, 60 s; units 2 and 3 partly follow unit 1
            rng(9); dur = 60;
            st = []; su = [];
            base = sort(rand(600, 1) * dur);
            for u = 1:n
                t = sort(rand(poissrnd((4 + 8 * rand) * dur), 1) * dur);
                if u <= 3, t = sort([t; base + 0.002 * u]); end
                st = [st; t]; su = [su; u * ones(numel(t), 1)]; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
        end

        function [W, T, B] = denseReference(spike_times, spike_units, n_units, duration, p)
            % The former dense implementation of computeConnectivityDDC.
            n_bins = ceil(duration / p.BinSize);
            V = zeros(n_bins, n_units);
            b = max(1, min(n_bins, round(spike_times / p.BinSize) + 1));
            V(sub2ind(size(V), b, spike_units)) = 1;
            W = test_ConnectivityDDC.denseCore(V, p.Threshold);
            S = zeros(n_units, n_units, p.N_Surrogates);
            for s = 1:p.N_Surrogates
                Vs = V;
                for u = 1:n_units
                    Vs(:, u) = circshift(V(:, u), randi(n_bins));
                end
                S(:, :, s) = test_ConnectivityDDC.denseCore(Vs, p.Threshold);
            end
            T = prctile(abs(S), p.Percentile, 3);
            B = double(abs(W) > T);
            B(1:n_units+1:end) = 0;
        end

        function W = denseCore(V, thr)
            N = size(V, 2);
            Fx = V - thr; Fx(Fx < 0) = 0;
            tmp = cov([Fx V]); Bm = tmp(1:N, N+1:end);
            dV = (-1/2 * V(1:end-2, :) + 1/2 * V(3:end, :));
            dV = [mean(dV); dV; mean(dV)];
            tmp = cov([dV V]); dCov = tmp(1:N, N+1:end);
            W = dCov / Bm;
        end

    end
end
