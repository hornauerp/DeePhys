classdef test_IntegrityFixes < matlab.unittest.TestCase
% TEST_INTEGRITYFIXES  Culture keys, responder correction, STTC surrogates, label propagation,
% sorting helpers, unit alignment, feature groups, normalisation, pooling and loading.
%
% Run: runtests('test_IntegrityFixes')

    properties
        Root
    end

    methods (TestMethodSetup)
        function makeRoot(tc)
            tc.Root = string(tempname);
            mkdir(tc.Root);
            tc.addTeardown(@() rmdir(tc.Root, 's'));
            w = warning('off', 'all');
            tc.addTeardown(@() warning(w));
        end
    end

    methods (Test)

        function testDatetimePlatingDateKeepsCultureKey(tc)
            sp = make_mock_sorting(tc.Root, 6, 300, 20000, 120);
            md = struct('ChipID', 'C1', 'PlatingDate', datetime(2024, 1, 15));
            sd = SpikeData.fromKilosort(sp, md);
            p = struct(); p.QC.N_Units = 1; p.QC.RPV = 1; p.Analyses.Catch22 = false;
            proc = RecordingProcessor(sd, p);
            proc.runQC(); proc.computeUnitFeatures();
            fs = FeatureStore.fromProcessors(proc);
            tc.verifyEqual(fs.MetadataTable.PlatingDate, "2024-01-15");
            keys = ["ChipID", "PlatingDate"];
            % the FeatureStore and the processor's own metadata give the same culture ID
            id_fs   = FeatureStore.buildCultureIDs(fs.MetadataTable, keys);
            id_proc = FeatureStore.buildCultureIDs(struct2table(sd.Metadata, 'AsArray', true), keys);
            tc.verifyEqual(id_fs, "C1_2024-01-15");
            tc.verifyEqual(id_proc, id_fs);
        end

        function testMissingCultureKeyWarns(tc)
            warning('on', 'FeatureStore:missingCultureKey');
            meta = table("C1", 'VariableNames', "ChipID");
            tc.verifyWarning(@() FeatureStore.buildCultureIDs(meta, ["ChipID", "PlatingDate"]), ...
                'FeatureStore:missingCultureKey');
        end

        function testNetworkCorrectionIsRateIndependent(tc)
            % every unit's rate x1.6 (a network gain) plus one unit x5: only that unit responds
            rng(3);
            n = 200; n_bins = 60;
            rate = exp(randn(n, 1)) * 2;                 % lognormal counts per bin
            pre  = poissrnd(repmat(rate, 1, n_bins));
            gain = 1.6 * ones(n, 1); gain(1) = 8;
            post = poissrnd(repmat(rate .* gain, 1, n_bins));
            r = bootstrapFiringRateResponse(pre, post, 500, 1e-4, 'NetworkCorrection', true);
            tc.verifyTrue(ismember(1, r.increase));
            % a common subtraction in counts called most fast units "increase" and most slow
            % units "decrease"
            tc.verifyLessThanOrEqual(numel(setdiff(r.increase, 1)) + numel(r.decrease), 0.05 * n);
        end

        function testSTTCJitterFollowsMaxLag(tc)
            p = RecordingProcessor.returnDefaultParams();
            tc.verifyEmpty(p.STTC.JitterSD);
            [cells, st, su] = test_IntegrityFixes.poissonUnits([1 2 3], 3, 60);
            params = struct('MaxLag', 0.1, 'SurrogateMethod', 'jitter', 'JitterSD', 0.05, ...
                'N_Surrogates', 2, 'Percentile', 95);
            warning('on', 'computeConnectivitySTTC:jitterTooSmall');
            tc.verifyWarning(@() computeConnectivitySTTC(cells, st, su, 3, 60, 20000, params), ...
                'computeConnectivitySTTC:jitterTooSmall');
        end

        function testSTTCSurrogatesSkipSilentUnits(tc)
            % units 2 and 4 are silent: every other unit must get its own surrogate
            [cells, st, su] = test_IntegrityFixes.poissonUnits([1 3 5], 5, 120);
            params = struct('MaxLag', 0.05, 'SurrogateMethod', 'jitter', 'JitterSD', [], 'EdgeTest', 'percentile', ...
                'N_Surrogates', 5, 'Percentile', 95);
            conn = computeConnectivitySTTC(cells, st, su, 5, 120, 20000, params);
            tc.verifyTrue(all(isfinite(conn.threshold([1 3], 5))));
            tc.verifyTrue(all(isnan(conn.threshold([1 3 5], 4))));
        end

        function testMissingLabelsAreNotCounterexamples(tc)
            ctc = CellTypeClassifier();
            ctc.Parameters = CellTypeClassifier.returnDefaultParams();
            ctc.Parameters.Bootstrap.GroundTruthMethod    = "metadata";
            ctc.Parameters.Bootstrap.LabelField           = "CellType";
            ctc.Parameters.Bootstrap.ResponsiveClassValue = "I";
            fs = FeatureStore();
            fs.UnitTable = table(categorical(["I"; "E"; missing; "E"]), 'VariableNames', "CellType");
            ctc.FeatureStore = fs;
            ctc.identifyResponsiveUnits();
            tc.verifyEqual(ctc.ResponsiveUnitIdx, [true false false false]);
            tc.verifyEqual(ctc.CounterexampleUnitIdx, [false true false true]);
        end

        function testPropagationLeavesUnreachedUnlabelled(tc)
            % nodes 1-4 a chain with seeds at both ends; node 5 isolated; nodes 6-7 a pair
            % without seeds
            G = zeros(7);
            G(1,2) = 1; G(2,3) = 1; G(3,4) = 1; G(6,7) = 1;
            G = G + G';
            [F, unreached] = CellTypeClassifier.propagateLabels(G, [1 4], [1 2], 1000, 1e-10);
            tc.verifyEqual(unreached', [false false false false true true true]);
            tc.verifyEqual(F(2, :), [2/3 1/3], 'AbsTol', 1e-6);
        end

        function testPropagationBalancesSeedCounts(tc)
            % a star around node 1 with three excitatory and one inhibitory seed leaves node 1
            % undecided: a node does not turn excitatory because that class has more seeds
            G = zeros(5); G(1, 2:5) = 1; G = G + G';
            F = CellTypeClassifier.propagateLabels(G, 2:5, [1 1 1 2], 1000, 1e-10);
            tc.verifyEqual(F(1, :), [0.5 0.5], 'AbsTol', 1e-9);
        end

        function testSplitSingleWindowFollowsTimes(tc)
            sp = make_mock_sorting(tc.Root, 6, 300, 20000, 120);
            fid = fopen(fullfile(sp, 'params.py'), 'a'); fprintf(fid, 'n_samples = %d\n', 20000 * 120); fclose(fid);
            out1 = split_sortings(string(sp), [0 0.5]);
            n1 = numel(readNPY(fullfile(out1, "spike_times.npy")));
            out2 = split_sortings(string(sp), [0.5 2]);     % same folder name, other window
            tc.verifyEqual(out2, out1);
            st2 = readNPY(fullfile(out2, "spike_times.npy"));
            tc.verifyNotEqual(numel(st2), n1);
            sd = SpikeData.fromKilosort(out2, struct('ChipID', 'c', 'PlatingDate', '2024-01-01'));
            tc.verifyEqual(sd.Duration, 90, 'AbsTol', 1e-9);
        end

        function testSplitRewritesFolderWithoutMarker(tc)
            sp = make_mock_sorting(tc.Root, 6, 300, 20000, 120);
            fid = fopen(fullfile(sp, 'params.py'), 'a'); fprintf(fid, 'n_samples = %d\n', 20000 * 120); fclose(fid);
            out = split_sortings(string(sp), [0 0.5; 0.5 1]);
            delete(fullfile(out(2), "split_window.txt"));   % as written by an older version
            copyfile(fullfile(sp, "params.py"), fullfile(out(2), "params.py"));
            out = split_sortings(string(sp), [0 0.5; 0.5 1]);
            tc.verifyTrue(isfile(fullfile(out(2), "split_window.txt")));
            sd = SpikeData.fromKilosort(out(2), struct('ChipID', 'c', 'PlatingDate', '2024-01-01'));
            tc.verifyEqual(sd.Duration, 30, 'AbsTol', 1e-9);
        end

        function testSortingPathListTrailingStar(tc)
            r = fullfile(tc.Root, 'tree');
            for w = 0:2, mkdir(fullfile(r, 'T001', 'Network', sprintf('w%03d', w))); end
            mkdir(fullfile(r, 'T001', 'Network', 'x000'));
            a = generate_sorting_path_list(string(r), {'T00*', 'Network', '*'});
            tc.verifyEqual(numel(a), 4);
            tc.verifyTrue(all(endsWith(a, ["w000", "w001", "w002", "x000"])));
        end

        function testUnitDataOrderIsChecked(tc)
            [p1, p2] = test_DataLayerFixes.twoProcessors(tc.Root);
            fs = FeatureStore.fromProcessors([p1 p2]);
            ud = [p2.Units, p1.Units];                      % same count, other order
            tc.verifyError(@() CellTypeClassifier(fs, ud), 'CellTypeClassifier:unitMisaligned');
        end

        function testActivityGroupSelectsAllActivityFeatures(tc)
            p1 = test_DataLayerFixes.processor(tc.Root, "ISIN");
            fs = FeatureStore.fromProcessors(p1);
            X = fs.unitMatrix("ActivityFeatures");
            tc.verifyTrue(ismember("ISIBimodality", string(X.Properties.VariableNames)));
            tc.verifyEqual(sort(string(X.Properties.VariableNames)), sort(unitFeatureNames("act")));
        end

        function testConstantFeatureInGroupIsZero(tc)
            rng(1);
            X = [randn(10, 1), [5 * ones(5, 1); randn(5, 1)]];   % column 2 constant in group a
            g = [repmat("a", 5, 1); repmat("b", 5, 1)];
            np = NormalizationPipeline.groupThenGlobal();
            [X_out, np_fit] = np.fit_transform(X, g);
            tc.verifyFalse(any(isnan(X_out(:))));
            tc.verifyLessThan(max(abs(X_out(1:5, 2))), 1e-12);
            X_new = np_fit.transform(X, g);
            tc.verifyEqual(X_new, X_out, 'AbsTol', 1e-12);
        end

        function testClipKeepsMissingValues(tc)
            X = [randn(9, 2); NaN 1];
            np = NormalizationPipeline.globalOnly();
            X_out = np.fit_transform(X, ones(10, 1));
            tc.verifyTrue(isnan(X_out(10, 1)));
        end

        function testPoolingUsesExactValues(tc)
            [idx, lbl] = MLPipeline.poolMetadataValues([1 2 3 1]', ["HET", "HET_KD", "WT"], "HET");
            tc.verifyEqual(idx', [1 2 2 1]);
            tc.verifyEqual(lbl, ["HET", "HET_KD/WT"]);
        end

        function testLoadFillsNewParameters(tc)
            proc = test_DataLayerFixes.processor(tc.Root, "ISIN");
            proc.Parameters.QC.N_Units = 3;                  % a stored, non-default value
            proc.Parameters = rmfield(proc.Parameters, 'STTC');      % as if added later
            proc.Parameters.QC = rmfield(proc.Parameters.QC, 'MinSpikeCount');
            f = fullfile(tc.Root, "proc.mat");
            proc.save(f);
            loaded = RecordingProcessor.load(f);
            defaults = RecordingProcessor.returnDefaultParams();
            tc.verifyEqual(loaded.Parameters.STTC, defaults.STTC);
            tc.verifyEqual(loaded.Parameters.QC.MinSpikeCount, defaults.QC.MinSpikeCount);
            tc.verifyEqual(loaded.Parameters.QC.N_Units, 3);
        end

    end

    methods (Static)

        function [cells, st, su] = poissonUnits(active, n_units, duration)
        % Independent 5 Hz Poisson trains for the active units; the others are silent.
            cells = repmat({zeros(0, 1)}, 1, n_units);
            st = []; su = [];
            for u = active
                t = sort(rand(round(5 * duration), 1) * duration);
                cells{u} = t;
                st = [st; t]; su = [su; u * ones(numel(t), 1)]; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
        end

    end
end
