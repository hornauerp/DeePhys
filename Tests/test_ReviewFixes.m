classdef test_ReviewFixes < matlab.unittest.TestCase
% TEST_REVIEWFIXES  Regression tests for the fourth review: pure-culture enrichment, segment
% QC, unit-table re-runs, drift-corrected responder tests, graph / ML / classifier state.
%
% Run: runtests('test_ReviewFixes')

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

        % ---- cell-type classifier --------------------------------------------------------

        function testPureCultureEnrichmentIsPooled(tc)
            % culture A: all labelled inhibitory, culture B: all labelled excitatory. Within a
            % culture the responder count is fixed, so a stratified test finds nothing.
            rng(5);
            X = [randn(100, 2) + [4 0]; randn(100, 2) - [4 0]];
            G = CellTypeClassifier.knnGraph(X, 10);
            resp = [true(100, 1); false(100, 1)];
            elig = true(200, 1);
            cult = [repmat("A", 100, 1); repmat("B", 100, 1)];
            p_comm = struct('LouvainResolution', 1, 'LouvainRestarts', 10, 'SeedFDR', 0.05, ...
                'Algorithm', "leiden", 'StratifyEnrichment', true);
            [inh, exc] = CellTypeClassifier.enrichmentSeeds(G, resp, elig, p_comm, 1, cult);
            tc.verifyNotEmpty(inh);
            tc.verifyTrue(all(inh <= 100));
            tc.verifyTrue(all(exc > 100));
        end

        function testCompositionCountsKnownLabels(tc)
            % a group made only of known seeds is their label fraction; known units are not
            % part of the EM of the other units
            rng(2);
            n = 400; y = rand(n, 1) < 0.2;
            Z = randn(n, 1) + 2.5 * y;
            train = [find(y, 40); find(~y, 40)];
            groups = repmat("g", n, 1); groups(train(1:10)) = "seedonly"; groups(train(41:50)) = "seedonly";
            known = nan(n, 1); known(train) = 1 + y(train);
            T = estimateClassComposition(Z, train, 1 + y(train), groups, NBoot=20, KnownLabels=known);
            tc.verifyEqual(T.Fraction(T.Group == "seedonly"), 0.5, 'AbsTol', 1e-12);
            g = T.Group == "g";
            truth = mean(y(groups == "g"));
            tc.verifyEqual(T.Fraction(g), truth, 'AbsTol', 0.05);
        end

        function testClearCacheDropsComposition(tc)
            ctc = CellTypeClassifier();
            ctc.Composition = struct('Probability', [0.2 0.9]);
            ctc.clearCache();
            tc.verifyFalse(isfield(ctc.Composition, 'Probability'));
        end

        function testStaleEmbeddingIsNotReused(tc)
            ctc = CellTypeClassifier();
            ctc.Parameters = CellTypeClassifier.returnDefaultParams();
            ctc.Parameters.Graph.Method = "exact";
            p = ctc.Parameters.UMAP; p.ComputeEmbedding = false; p.AutoNNeighbors = false; p.NNeighbors = 5;
            X = randn(30, 4);
            ctc.UMAP = struct('stale', true);
            ctc.Reduction = struct('Unsupervised', zeros(30, 2));    % same N, other features
            [~, red] = evalc('ctc.runUnsupervisedUMAP(X, p, true)');
            tc.verifyEmpty(red);
        end

        function testMissingFeaturesIgnoredByGroupZScore(tc)
            X = [1 10; 2 20; 3 30; NaN NaN];
            Z = zscoreByGroup(X, ones(4, 1), true(4, 1), 1);
            tc.verifyEqual(Z(1:3, 1), [-1; 0; 1], 'AbsTol', 1e-12);
            tc.verifyTrue(all(isnan(Z(4, :))));
        end

        function testQualityReportUsesCallDirection(tc)
            ctc = test_ReviewFixes.reportClassifier();
            ctc.Parameters.Bootstrap.Direction = "decrease";
            ctc.ResponsivenessDetail = struct('calibration', struct('fdr_inc', 0.9, 'fdr_dec', 0.05, ...
                'rate_inc_treated', 0.01, 'rate_inc_control', 0.01, ...
                'rate_dec_treated', 0.2, 'rate_dec_control', 0.01));
            [~, r] = evalc('ctc.qualityReport()');
            row = r.checks(r.checks.Check == "Responder FDR", :);
            tc.verifyEqual(row.Status, "PASS");
            tc.verifyTrue(startsWith(row.Detail, "decrease"));
            % with response rules, every rule is checked in its own direction
            ctc.ResponsivenessDetail = struct('rules', struct('Direction', {"increase", "decrease"}, ...
                'Calibration', {struct('fdr_inc', 0.05, 'fdr_dec', 0.9, 'rate_inc_treated', 0.2, ...
                'rate_inc_control', 0.01, 'rate_dec_treated', 0.01, 'rate_dec_control', 0.01), ...
                struct('fdr_inc', 0.05, 'fdr_dec', 0.6, 'rate_inc_treated', 0.2, ...
                'rate_inc_control', 0.01, 'rate_dec_treated', 0.02, 'rate_dec_control', 0.012)}));
            [~, r] = evalc('ctc.qualityReport()');
            tc.verifyEqual(r.checks.Status(r.checks.Check == "Responder FDR"), "WARN");
        end

        function testCalibratedPValueIsDriftCorrected(tc)
            % treated culture A: every unit at half its rate after the drug, no responders
            [ud, meta] = test_ReviewFixes.driftExperiment(0.5);
            p = struct('BinSize', 20, 'HighDoses', 2);
            T = CellTypeClassifier.responseTable(ud, meta, ["ChipID", "PlatingDate"], ...
                "Concentration", 0, ["A_1"; "B_1"], p);
            is_a = startsWith(T.UnitID, "A");
            tc.verifyLessThan(mean(T.p(is_a) < 0.01), 0.1);
        end

        function testDiscardSegmentUnitsKeepsOriginal(tc)
            [ud, meta] = test_ReviewFixes.driftExperiment(1);
            fs = FeatureStore();
            fs.MetadataTable = meta;
            fs.UnitTable = table(string({ud.UnitID})', string({ud.RecordingID})', ...
                'VariableNames', {'UnitID', 'RecordingID'});
            ctc = CellTypeClassifier();
            ctc.Parameters = CellTypeClassifier.returnDefaultParams();
            ctc.Parameters.UMAP.GroupingVar = "Concentration";
            ctc.Parameters.UMAP.GroupingValues = 0;
            ctc.FeatureStore = fs; ctc.UnitDataArray = ud;
            evalc('ctc.discardSegmentUnits()');
            n_full = numel(ctc.OriginalUnitDataArray);
            evalc('ctc.discardSegmentUnits()');                % a re-run cell
            tc.verifyEqual(n_full, numel(ud));
            tc.verifyEqual(numel(ctc.OriginalUnitDataArray), numel(ud));
        end

        % ---- data layer ------------------------------------------------------------------

        function testSegmentQCUsesParentSorting(tc)
            % template 0 fires only in the first half: the second window keeps it (QC on the
            % parent sorting), so every segment holds the same units
            sp = test_ReviewFixes.halfSilentSorting(tc.Root);
            out = split_sortings(string(sp), [0 5; 5 10]);
            md = struct('ChipID', 'c', 'PlatingDate', '2024-01-01');
            p.QC.N_Units = 1; p.Analyses.Regularity = false; p.Analyses.Catch22 = false;
            ids = cell(1, 2);
            for w = 1:2
                proc = RecordingProcessor(SpikeData.fromKilosort(out(w), md), p);
                evalc('proc.runQC()');
                ids{w} = [proc.Units.TemplateID];
            end
            tc.verifyEqual(ids{2}, ids{1});
            tc.verifyTrue(ismember(1, ids{2}));
            evalc('proc.computeUnitFeatures()');            % a unit without spikes in the window
            tc.verifyEqual(height(proc.UnitFeatureTable), numel(ids{2}));
            % window-only QC drops the silent unit
            p.QC.RateReference = "recording";
            proc = RecordingProcessor(SpikeData.fromKilosort(out(2), md), p);
            evalc('proc.runQC()');
            tc.verifyFalse(ismember(1, [proc.Units.TemplateID]));
        end

        function testRecomputeUnitFeaturesKeepsOtherSteps(tc)
            sp = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            p.QC.N_Units = 1; p.Analyses.Regularity = false; p.Analyses.Catch22 = false;
            proc = RecordingProcessor(SpikeData.fromKilosort(sp, struct('ChipID', 'c', 'PlatingDate', '2024-01-01')), p);
            evalc('proc.runQC(); proc.computeUnitFeatures()');
            n = height(proc.UnitFeatureTable);
            tc.assertGreaterThan(n, 0);
            proc.UnitFeatureTable.Parent_ACG1 = (1:n)';
            proc.UnitFeatureTable.BurstFractionSpikes = (1:n)' / 10;
            proc.Status.UnitFeatures = "pending";
            evalc('proc.computeUnitFeatures("ActivityFeatures")');
            cols = string(proc.UnitFeatureTable.Properties.VariableNames);
            tc.verifyEqual(proc.UnitFeatureTable.Parent_ACG1, (1:n)');
            tc.verifyEqual(proc.UnitFeatureTable.BurstFractionSpikes, (1:n)' / 10);
            tc.verifyFalse(any(cols == "HalfWidth"));        % the step's own stale columns go
            tc.verifyFalse(any(cols == "ACG1"));
            tc.verifyTrue(any(cols == "FiringRate"));
        end

        function testGoodUnitsAcceptsIndices(tc)
            sp = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            sd = SpikeData.fromKilosort(sp, struct('ChipID', 'c', 'PlatingDate', '2024-01-01'));
            qc = @(gu) RecordingProcessor(sd, struct('QC', struct('GoodUnits', gu, 'N_Units', 1)));
            proc = qc([3 5]);                                   % template numbers
            evalc('proc.runQC()');
            tc.verifyEqual([proc.Units.TemplateID], [3 5]);
            proc = qc(logical([1 0 1 0 0]));                    % mask
            evalc('proc.runQC()');
            tc.verifyEqual([proc.Units.TemplateID], [1 3]);
            proc = qc(logical([1 0]));                          % wrong length
            tc.verifyError(@() proc.runQC(), 'RecordingProcessor:goodUnits');
            proc = qc(7);                                       % no such template
            tc.verifyError(@() proc.runQC(), 'RecordingProcessor:goodUnits');
        end

        function testSilentWindowsGetDistinctIDs(tc)
            a = SpikeData.contentID([], [], 20000, "/x/sorting/w1");
            b = SpikeData.contentID([], [], 20000, "/x/sorting/w2");
            tc.verifyNotEqual(a, b);
            % recordings with spikes stay path-independent
            tc.verifyEqual(SpikeData.contentID(1:5, ones(1, 5), 20000, "/a"), ...
                SpikeData.contentID(1:5, ones(1, 5), 20000, "/b"));
        end

        function testRegistryReadsTextDIV(tc)
            tc.assumeTrue(RecordingDatabase.isAvailable(), 'Database Toolbox not available.');
            sp = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            sd = SpikeData.fromKilosort(sp, struct('ChipID', 'DIVChip', 'PlatingDate', '240101', 'DIV', '14'));
            p.Analyses.Catch22 = false; p.QC.N_Units = 1;
            proc = RecordingProcessor(sd, p);
            evalc('proc.runQC()');
            f = fullfile(tc.Root, 'reg.db');
            RecordingDatabase.resetInstance();
            db = RecordingDatabase.instance(f);
            tc.addTeardown(@() RecordingDatabase.resetInstance());
            db.registerRecording(proc);
            T = db.queryRecordings(ChipID="DIVChip");
            tc.verifyEqual(double(T.div(1)), 14);
        end

        function testMetadataDatetimeIsISO(tc)
            % an Excel lookup with a date column is read as datetime
            xlsx = fullfile(tc.Root, 'lookup.xlsx');
            writetable(table("P001_03", datetime(2024, 1, 1), "WT", ...
                'VariableNames', {'Chip_IDs', 'PlatingDate', 'Mutation'}), xlsx);
            T = generateMetadataTable("/r/240115/P001/Network/well003", xlsx, [3 4 6]);
            tc.verifyEqual(string(T.PlatingDate(1)), "2024-01-01");
            tc.verifyEqual(T.DIV(1), 14);
        end

        % ---- features --------------------------------------------------------------------

        function testSmallWorldnessWithIsolatedNode(tc)
            rng(3);
            A = double(rand(40) < 0.15); A = triu(A, 1); A = A + A';
            A(40, :) = 0; A(:, 40) = 0;                          % one isolated unit
            conn.STTC.bu = A;
            nw = computeGraphFeatures(conn, struct('SmallWorldnessIterations', 5));
            sw = nw{1, startsWith(string(nw.Properties.VariableNames), "SmallWorldness")};
            tc.verifyTrue(all(isfinite(sw) & sw > 0));
        end

        function testCatch22HasNoRateFeatures(tc)
            tc.assumeTrue(exist('catch22_DN_HistogramMode_5', 'file') == 3, 'catch22 not compiled');
            T = computeCatch22(sort(rand(500, 1) * 60), 60, 0.1);
            tc.verifyEqual(width(T), 22);
            tc.verifyFalse(any(ismember(T.Properties.VariableNames, {'DN_Mean', 'DN_Spread_Std'})));
        end

        % ---- ML --------------------------------------------------------------------------

        function testCombatWithoutCovariates(tc)
            rng(1);
            X = randn(60, 4) + repelem([0; 3], 30);
            g = repelem([1; 2], 30);
            train = true(60, 1); train(1:5:end) = false;
            [Xtr, Xte] = MLUtils.normalizeFeatures(X, train, ~train, g, "combat", []);
            tc.verifySize(Xtr, [sum(train) 4]);
            tc.verifySize(Xte, [sum(~train) 4]);
            tc.verifyError(@() MLUtils.normalizeFeatures(X, train, ~train, [], "combat", []), ...
                'MLUtils:combatNeedsGroups');
        end

        function testRandomForestSeedIsReproducible(tc)
            rng(1);
            X = randn(120, 6); Y = double(X(:, 1) + 0.5 * randn(120, 1) > 0) + 1;
            grp = repelem((1:12)', 10);
            opts = struct('Algorithm', "rf", 'KFold', 3, 'Seed', 5, 'CVGroups', grp);
            [~, r1] = evalc('Classifier.classify(X, Y, opts)');
            [~, r2] = evalc('Classifier.classify(X, Y, opts)');
            tc.verifyEqual(vertcat(r1.scores), vertcat(r2.scores));
        end

        function testValidationReportsBothClasses(tc)
            pred = [ones(1, 90), 2 * ones(1, 10), ones(1, 5), 2 * ones(1, 5)];
            truth = [ones(1, 100), 2 * ones(1, 10)];
            [~, m] = evalc('validateClassification(pred, truth)');
            close all
            tc.verifyEqual(m.PerClassRecall, [0.9 0.5], 'AbsTol', 1e-12);
            tc.verifyEqual(m.Recall, 0.7, 'AbsTol', 1e-12);
        end

        % ---- figure scripts --------------------------------------------------------------

        function testFigureScriptsFindTheirData(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            files = dir(fullfile(root, 'Figures', 'Paper', '*.m'));
            for f = files'
                txt = fileread(fullfile(f.folder, f.name));
                tok = regexp(txt, 'fullfile\(deephys_root, ("[^)]*")\)', 'tokens');
                for k = 1:numel(tok)
                    parts = eval(['[' tok{k}{1} ']']);
                    if any(parts == "Plots"), continue, end
                    d = fullfile(root, parts{:});
                    tc.verifyTrue(isfolder(d), sprintf('%s: %s', f.name, d));
                end
            end
        end

    end

    methods (Static)

        function sp = halfSilentSorting(root)
            % 12 templates over 600 s; template 0 fires only in the first 300 s
            fs = 20000; n_units = 12; dur = 600;
            sp = make_mock_sorting(root, n_units, 10, fs, dur);
            rng(11);
            t = []; u = [];
            for k = 0:n_units - 1
                if k == 0
                    tk = rand(300, 1) * 290 + 5;
                else
                    tk = rand(300, 1) * (dur - 10) + 5;
                end
                t = [t; round(tk * fs)]; u = [u; repmat(k, 300, 1)]; %#ok<AGROW>
            end
            [t, o] = sort(t); u = u(o);
            writeNPY(uint64(t), fullfile(sp, 'spike_times.npy'));
            writeNPY(uint32(u), fullfile(sp, 'spike_templates.npy'));
            fid = fopen(fullfile(sp, 'params.py'), 'a'); fprintf(fid, 'n_samples = %d\n', fs * dur); fclose(fid);
        end

        function [ud, meta] = driftExperiment(gain)
            % Cultures A (treated) and B (control), recordings at doses 0, 1, 10 (600 s).
            % Every unit of A fires at gain x its rate at doses > 0; no single responders.
            rng(8);
            chips = ["A", "B"]; conc = [0 1 10]; dur = 600;
            RecordingID = strings(0, 1); ChipID = RecordingID; PlatingDate = RecordingID; Concentration = [];
            ud = UnitData.empty(1, 0);
            for c = chips
                base_rate = 2 + 3 * rand(50, 1);
                for d = conc
                    rid = c + "_" + d;
                    RecordingID(end+1, 1) = rid; ChipID(end+1, 1) = c; PlatingDate(end+1, 1) = "1"; %#ok<AGROW>
                    Concentration(end+1, 1) = d; %#ok<AGROW>
                    for u = 1:50
                        r = base_rate(u);
                        if c == "A" && d > 0, r = gain * r; end
                        st = sort(rand(poissrnd(r * dur), 1) * dur);
                        x = UnitData(u, u, st, zeros(10, 1), rid, 10000, dur);
                        x.UnitID = c + "_u" + u;
                        ud(end+1) = x; %#ok<AGROW>
                    end
                end
            end
            meta = table(RecordingID, ChipID, PlatingDate, Concentration);
        end

        function ctc = reportClassifier()
            % the minimum qualityReport reads: seeds, labels, units and metadata
            ctc = CellTypeClassifier();
            ctc.Parameters = CellTypeClassifier.returnDefaultParams();
            ctc.TrainLabels = struct('sorted_y_train', [ones(1, 40), 2 * ones(1, 40)]);
            ctc.UnitLabels = [ones(1, 8), 2 * ones(1, 2)];
            fs = FeatureStore();
            fs.UnitTable = table(repmat("r1", 10, 1), "u" + (1:10)', 'VariableNames', {'RecordingID', 'UnitID'});
            fs.MetadataTable = table("r1", "c1", "2024-01-01", 'VariableNames', {'RecordingID', 'ChipID', 'PlatingDate'});
            ctc.FeatureStore = fs;
        end

    end
end
