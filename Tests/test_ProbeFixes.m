classdef test_ProbeFixes < matlab.unittest.TestCase
% TEST_PROBEFIXES  Multiple-testing corrections with missing p-values, stratified grouped CV
% with text labels, ComBat with missing values, side-effect-free legacy loading,
% sub-sample waveform timing, burst statistics,
% log-ISI recruitment, burst propagation speed and the spatial FC decay fit.
%
% Run: runtests('test_ProbeFixes')

    methods (Test)

        function testFdrHolmIgnoreNaN(tc)
            % NaN p-values stay NaN and do not count as tests
            tc.verifyEqual(StatisticalTest.fdrCorrection([0.01 NaN 0.03]), [0.02; NaN; 0.03], 'AbsTol', 1e-12);
            tc.verifyEqual(StatisticalTest.holmCorrection([0.01 NaN 0.03]), [0.02; NaN; 0.03], 'AbsTol', 1e-12);
            % unchanged without NaN (BH step-up, Holm step-down)
            tc.verifyEqual(StatisticalTest.fdrCorrection([0.01 0.04 0.03 0.5]), [0.04; 0.04*4/3; 0.04*4/3; 0.5], 'AbsTol', 1e-12);
            tc.verifyEqual(StatisticalTest.holmCorrection([0.01 0.04 0.03 0.5]), [0.04; 0.09; 0.09; 0.5], 'AbsTol', 1e-12);
            tc.verifyEmpty(StatisticalTest.fdrCorrection([]));
        end

        function testGroupedCVStringLabels(tc)
            % stratification with text labels (metadata such as Mutation is often a string)
            g = repelem("rec" + (1:10)', 5);
            y = repelem(["WT"; "KO"], 25);
            gcv = GroupedCV.byGroups(g, 5, y);
            tc.verifyEqual(gcv.NumFolds, 5);
            gcv = GroupedCV.byGroups(g, 5, categorical(y));
            tc.verifyEqual(gcv.NumFolds, 5);
            gcv = GroupedCV.byGroups(g, 5, cellstr(y));
            tc.verifyEqual(gcv.NumFolds, 5);
        end

        function testComBatKeepsMissingLocal(tc)
            rng(1);
            b = repelem(["a"; "b"; "c"], 12);
            X = randn(36, 4) + (b == "b") * 2;
            np = NormalizationPipeline.combatThenGlobal();
            Y0 = np.fit_transform(X, b);
            X1 = X; X1(5, 2) = NaN;                       % one missing entry
            Y1 = np.fit_transform(X1, b);
            tc.verifyEqual(sum(isnan(Y1(:))), 1);
            tc.verifyTrue(isnan(Y1(5, 2)));
            tc.verifyLessThan(max(abs(Y1(:, [1 3 4]) - Y0(:, [1 3 4])), [], 'all'), 0.05);
            X2 = X; X2(:, 3) = NaN;                       % one all-missing feature
            Y2 = np.fit_transform(X2, b);
            tc.verifyEqual(find(all(isnan(Y2), 1)), 3);
            tc.verifyFalse(any(isnan(Y2(:, [1 2 4])), 'all'));
        end

        function testLegacyLoadWritesNoDatabase(tc)
            db_file = fullfile(tempdir, "probe_fix_registry_" + string(feature('getpid')) + ".db");
            if isfile(db_file), delete(db_file); end
            old = getenv('DEEPHYS_DB_PATH');
            setenv('DEEPHYS_DB_PATH', db_file);
            RecordingDatabase.instance("", true);
            c = onCleanup(@() restoreDB(old, db_file));
            obj = MEArecording(); %#ok<NASGU>
            f = [tempname '.mat'];
            save(f, 'obj');
            S = load(f, 'obj');
            delete(f);
            tc.verifyClass(S.obj, 'MEArecording');
            tc.verifyFalse(isfile(db_file), 'loading a legacy recording must not create a registry');
        end

        function testWaveformTimingBetweenSamples(tc)
            % 10 kHz: trough-to-peak delays between samples are no longer rounded to 0.1 ms
            sr = 1e4; t = (0:81)' / sr * 1e3;
            D = [0.32 0.36 0.42 0.46];
            est = zeros(size(D));
            for i = 1:numel(D)
                w = -exp(-((t - 2) / 0.12).^2) + 0.35 * exp(-((t - 2 - D(i)) / 0.25).^2);
                f = computeWaveformFeatures(w / max(abs(w)), sr);
                est(i) = f{1}.T2Pdelay;
            end
            tc.verifyLessThan(max(abs(est - D)), 0.03);
            tc.verifyEqual(numel(unique(est)), numel(D));
        end

        function testBurstSpreadUsesFilteredBursts(tc)
            % 20 bursts of ~0.2 s and one 200 s outlier: the SD follows the kept bursts
            starts = (0:20)' * 20;
            ends = starts + 0.2 + 0.01 * (0:20)';
            ends(end) = starts(end) + 200;
            st = sort([starts + 0.05; starts + 0.1; ends - 0.01]);
            su = mod((1:numel(st))', 5) + 1;
            p = struct('Binning', 0.1, 'Outlier', struct('Method', "median", 'ThresholdFactor', 5));
            F = computeBurstStatistics(struct('T_start', starts', 'T_end', ends'), st, su, 5, p);
            tc.verifyLessThan(F.StdBurstDuration, 1);
            tc.verifyLessThan(F.BurstDuration, 1);
        end

        function testLogISIRequiresRecruitment(tc)
            % independent unit bursts and no network bursts: log-ISI on the pooled train finds
            % many events; recruitment removes them. Real network bursts are kept.
            w = warning('off', 'all'); c = onCleanup(@() warning(w));
            p = struct('Method', "logISI", 'Binning', 0.1, 'MergeFactor', 0.5, ...
                'Outlier', struct('Method', "", 'ThresholdFactor', 5));
            [st, su] = simulateNetworkBursts(Mode="independent", Seed=2);
            [fr, cv] = test_ProbeFixes.unitStats(st, su, 300);
            b_ind = computeBursts(st, su, fr, cv, numel(fr), 300, p);
            [st, su, truth] = simulateNetworkBursts(Mode="network", Seed=2);
            [fr, cv] = test_ProbeFixes.unitStats(st, su, 300);
            b_net = computeBursts(st, su, fr, cv, numel(fr), 300, p);
            tc.verifyLessThan(numel(b_ind.T_start), 20);     % was > 1000 without recruitment
            mid = (b_net.T_start(:) + b_net.T_end(:)) / 2;
            found = arrayfun(@(k) any(mid >= truth.T_start(k) - 0.1 & mid <= truth.T_end(k) + 0.1), ...
                1:numel(truth.T_start));
            tc.verifyGreaterThanOrEqual(mean(found), 0.9);    % true network bursts are kept
            tc.verifyLessThan(numel(b_net.T_start), 5 * numel(truth.T_start));
        end

        function testPropagationSpeedFromMedianSlope(tc)
            % a wave at 50 um/ms in 9 of 10 bursts, one burst with a flat slope
            n = 20; xy = [(1:n)' * 100, zeros(n, 1)];
            units = arrayfun(@(k) UnitData(k, k, [], []), 1:n);
            st = []; su = [];
            for b = 1:10
                t0 = (b - 1) * 10;
                if b == 10, lat = 0.005 + 1e-6 * (1:n)';      % near-flat: 1/slope explodes
                else,       lat = (0:n-1)' * 100 / 50 / 1000; end   % 50 um/ms
                st = [st; t0 + lat]; su = [su; (1:n)']; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
            bursts = struct('T_start', (0:9) * 10, 'T_end', (0:9) * 10 + 1);
            T = computeSpatialBurstFeatures(xy, units, bursts, st, su, [], struct());
            tc.verifyEqual(T.MeanBurstPropagationSpeed, 50, 'RelTol', 0.05);
            % latencies unrelated to distance: no speed
            rng(4);
            st = []; su = [];
            for b = 1:10
                st = [st; (b - 1) * 10 + 0.05 * rand(n, 1)]; su = [su; (1:n)']; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
            T0 = computeSpatialBurstFeatures(xy, units, bursts, st, su, [], struct());
            tc.verifyTrue(isnan(T0.MeanBurstPropagationSpeed));
        end

        function testSpatialDecayFit(tc)
            % STTC = 0.3 exp(-d / 300 um) - 0.02 + noise (negative values included)
            rng(3); n = 80;
            units = arrayfun(@(k) UnitData(k, k, [], []), 1:n);
            xy = rand(n, 2) .* [3800 2100];
            D = squareform(pdist(xy));
            W = 0.3 * exp(-D / 300) - 0.02 + 0.03 * randn(n);
            W = (W + W') / 2; W(1:n+1:end) = NaN;
            T = computeSpatialConnectivityFeatures(xy, units, struct('STTC', struct('wu', W)));
            tc.verifyEqual(T.SpatialDecayTau, 300, 'RelTol', 0.15);
            % coupling that rises with distance does not decay: NaN
            W1 = 1e-4 * D + 0.01 * randn(n); W1 = (W1 + W1') / 2; W1(1:n+1:end) = NaN;
            T1 = computeSpatialConnectivityFeatures(xy, units, struct('STTC', struct('wu', W1)));
            tc.verifyTrue(isnan(T1.SpatialDecayTau));
            % decay far shorter than the closest pair (no structure at the measured distances): NaN
            W2 = 0.3 * exp(-D / 0.5) + 0.01 * randn(n); W2 = (W2 + W2') / 2; W2(1:n+1:end) = NaN;
            T2 = computeSpatialConnectivityFeatures(xy, units, struct('STTC', struct('wu', W2)));
            tc.verifyTrue(isnan(T2.SpatialDecayTau));
        end

    end

    methods (Static)
        function [fr, cv] = unitStats(st, su, dur)
            n = max(su);
            fr = accumarray(su, 1, [n 1]) / dur;
            cv = accumarray(su, st, [n 1], @(t) std(diff(sort(t))) / max(mean(diff(sort(t))), eps));
        end
    end
end

function restoreDB(old, db_file)
    RecordingDatabase.instance("", true);
    setenv('DEEPHYS_DB_PATH', old);
    if isfile(db_file), delete(db_file); end
end
