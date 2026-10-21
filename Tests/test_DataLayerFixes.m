classdef test_DataLayerFixes < matlab.unittest.TestCase
% TEST_DATALAYERFIXES  Table assembly, caching, step status and sorting-path helpers.
%
% Run: runtests('test_DataLayerFixes')

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

        function testStackedParentACGKeepsBinOrder(tc)
            [p1, p2] = test_DataLayerFixes.twoProcessors(tc.Root);
            n = numel(p2.Units);
            p2.UnitFeatureTable = [p2.UnitFeatureTable, ...
                array2table(repmat(1:401, n, 1), 'VariableNames', "Parent_ACG" + (1:401))];
            fs = FeatureStore.fromProcessors([p1 p2]);   % first recording has no parent ACG
            cols = string(fs.UnitTable.Properties.VariableNames);
            pc = cols(startsWith(cols, "Parent_ACG"));
            tc.verifyEqual(pc(1:3), ["Parent_ACG1", "Parent_ACG2", "Parent_ACG3"]);
            X = fs.unitMatrix('ACG', 'ACG');
            tc.verifyEqual(X{end, 1:5}, 1:5);
            % rows without a parent ACG keep their own ACG (not NaN)
            tc.verifyFalse(any(isnan(X{1, :})));
        end

        function testEmptyCustomSetSelectsNothing(tc)
            [p1, p2] = test_DataLayerFixes.twoProcessors(tc.Root);
            fs = FeatureStore.fromProcessors([p1 p2]);
            FeatureCatalog.registerSet("emptytestset", string.empty);
            tc.addTeardown(@() FeatureCatalog.removeSet("emptytestset"));
            X = fs.unitMatrix('all', FeatureSet="emptytestset");
            tc.verifyEqual(width(X), 0);
        end

        function testExtractionCacheFollowsUnitOrder(tc)
            base = (0.05:0.1:60)';                  % no ISI below 100 ms
            dbl  = sort([base; base + 0.005]);      % doublets at 5 ms
            u1 = UnitData(1, 1, base, -ones(10, 1), "recA", 10000, 60); u1.UnitID = "cacheA";
            u2 = UnitData(2, 2, dbl,  -ones(10, 1), "recA", 10000, 60); u2.UnitID = "cacheB";
            c1 = CellTypeClassifier(); c1.Parameters = CellTypeClassifier.returnDefaultParams();
            [~, a1] = c1.getOrExtract([u1 u2]);
            c2 = CellTypeClassifier(); c2.Parameters = c1.Parameters;
            [~, a2] = c2.getOrExtract([u2 u1]);     % disk cache of the first call
            [~, a3] = c1.getOrExtract([u2 u1]);     % in-memory cache of the first call
            tc.verifyEqual(a2, a1(:, [2 1]));
            tc.verifyEqual(a3, a1(:, [2 1]));
        end

        function testFailedSubAnalysisMarksStepFailed(tc)
            proc = test_DataLayerFixes.processor(tc.Root, "no_such_method");
            proc.computeNetworkFeatures();
            tc.verifyEqual(proc.Status.NetworkFeatures, "failed");
            proc.Parameters.Bursts.Method = "ISIN";
            proc.computeNetworkFeatures();          % re-runs because the step failed
            tc.verifyEqual(proc.Status.NetworkFeatures, "done");
            tc.verifyTrue(any(string(proc.NetworkFeatureTable.Properties.VariableNames) == "BurstDuration"));
        end

        function testConnectivityAddsRequestedMethod(tc)
            proc = test_DataLayerFixes.processor(tc.Root, "ISIN");
            proc.computeNetworkFeatures();
            proc.computeConnectivity();             % CCG
            n_cols = width(proc.NetworkFeatureTable);
            proc.computeConnectivity("STTC");
            tc.verifyTrue(isfield(proc.Connectivity, 'STTC'));
            names = string(proc.NetworkFeatureTable.Properties.VariableNames);
            tc.verifyEqual(numel(unique(names)), numel(names));   % graph columns replaced, not duplicated
            tc.verifyGreaterThan(width(proc.NetworkFeatureTable), n_cols);
        end

        function testQCKeepsTemplateLength(tc)
            sp = make_mock_sorting(tc.Root, 6, 300, 20000, 120);
            sd = SpikeData.fromKilosort(sp, struct('ChipID', 'C1', 'PlatingDate', '2024-01-01'));
            n_samp = size(sd.TemplateWaveforms, 2);
            sd.TemplateWaveforms(:, 1, :) = 0;      % a sample that is zero in every template
            proc = RecordingProcessor(sd, struct('QC', struct('N_Units', 1, 'RPV', 1)));
            proc.runQC();
            tc.verifyEqual(numel(proc.Units(1).ReferenceWaveform), n_samp);
        end

        function testSortingPathListWildcardAndLiteral(tc)
            r = fullfile(tc.Root, 'tree');
            for w = 0:2, mkdir(fullfile(r, 'T001', 'Network', sprintf('w%03d', w), 'sorter_output')); end
            a = generate_sorting_path_list(string(r), {'T00*', 'Network', 'w*'});
            b = generate_sorting_path_list(string(r), {'T00*', 'Network', 'w*', 'sorter_output'});
            tc.verifyEqual(numel(a), 3);
            tc.verifyTrue(all(endsWith(a, ["w000", "w001", "w002"])));
            tc.verifyEqual(numel(b), 3);
            tc.verifyTrue(all(endsWith(b, "sorter_output")));
        end

        function testSplitWindowsDurationAndSortingID(tc)
            sp = make_mock_sorting(tc.Root, 6, 300, 20000, 120);
            fid = fopen(fullfile(sp, 'params.py'), 'a'); fprintf(fid, 'n_samples = %d\n', 20000 * 120); fclose(fid);
            out = split_sortings(string(sp), [0 0.5; 0.5 1.0]);
            md = struct('ChipID', 'c', 'PlatingDate', '2024-01-01');
            sd1 = SpikeData.fromKilosort(out(1), md);
            sd2 = SpikeData.fromKilosort(out(2), md);
            tc.verifyEqual(sd1.Duration, 30, 'AbsTol', 1e-9);
            tc.verifyEqual(sd1.SortingID, sd2.SortingID);
            tc.verifyNotEqual(sd1.SortingID, sd1.RecordingID);
        end

        function testAnalysisLogLeavesGlobalRNG(tc)
            rng(5); a = rand;
            rng(5); AnalysisLog(); b = rand;
            tc.verifyEqual(a, b);
        end

        function testTemplateBaseAmbiguityWarns(tc)
            d = fullfile(tc.Root, 'pbase'); mkdir(d);
            writeNPY(zeros(5, 4, 2, 'single'), fullfile(d, 'templates.npy'));
            warning('on', 'ParentSpikeLoader:indexBaseAmbiguous');
            tc.verifyWarning(@() ParentSpikeLoader.templateIndexBase(d, [1 2 3 4 1 2 3 4]'), ...
                'ParentSpikeLoader:indexBaseAmbiguous');
        end

        function testUnknownBurstMethodErrors(tc)
            tc.verifyError(@() computeBursts((0:0.1:10)', ones(101, 1), 1, 1, 1, 10, ...
                struct('Method', "bogus")), 'computeBursts:unknownMethod');
        end

    end

    methods (Static)

        function proc = processor(root, burst_method)
            sp = make_mock_sorting(root, 8, 300, 20000, 120);
            sd = SpikeData.fromKilosort(sp, struct('ChipID', 'C1', 'PlatingDate', '2024-01-01'));
            p = struct(); p.QC.N_Units = 1; p.QC.RPV = 1; p.Analyses.Catch22 = false;
            p.Analyses.Connectivity = "CCG"; p.GraphFeatures.SmallWorldnessIterations = 5;
            p.STTC.N_Surrogates = 5;
            p.Bursts.Method = burst_method;
            proc = RecordingProcessor(sd, p);
            proc.runQC(); proc.computeUnitFeatures();
        end

        function [p1, p2] = twoProcessors(root)
            p1 = test_DataLayerFixes.processor(root, "ISIN");
            p1.computeNetworkFeatures();
            p2 = test_DataLayerFixes.processor(root, "ISIN");   % own SpikeData object
            p2.computeNetworkFeatures();
            p2.SpikeData.RecordingID = "other_rec";
            for u = 1:numel(p2.Units), p2.Units(u).RecordingID = "other_rec"; end
        end

    end
end
