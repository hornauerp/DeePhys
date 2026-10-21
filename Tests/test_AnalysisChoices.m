classdef test_AnalysisChoices < matlab.unittest.TestCase
% TEST_ANALYSISCHOICES  STTC edge test, feature-definition versions, culture-stratified
% enrichment, the LvR constant and the burst-participation feature group.
%
% Run: runtests('test_AnalysisChoices')

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

        function testSTTCFdrEdgesOnIndependentUnits(tc)
            % independent units: the per-pair percentile leaves ~5% false edges, FDR almost none
            rng(2);
            n = 20; dur = 300;
            [cells, st, su] = test_AnalysisChoices.poissonTrains(n, dur, 2);
            p = struct('MaxLag', 0.05, 'SurrogateMethod', 'jitter', 'JitterSD', [], ...
                'N_Surrogates', 40, 'Percentile', 95, 'EdgeTest', 'fdr', 'FDR', 0.05);
            c_fdr = computeConnectivitySTTC(cells, st, su, n, dur, 20000, p);
            p.EdgeTest = 'percentile';
            c_pct = computeConnectivitySTTC(cells, st, su, n, dur, 20000, p);
            n_pairs = n * (n - 1) / 2;
            tc.verifyLessThanOrEqual(nnz(triu(c_fdr.bu, 1)), 1);
            tc.verifyGreaterThan(nnz(triu(c_pct.bu, 1)), 0.01 * n_pairs);
            tc.verifyEqual(c_fdr.bu, c_fdr.bu.');
            tc.verifyFalse(any(diag(c_fdr.bu)));
        end

        function testSTTCFdrKeepsSynchronousPair(tc)
            rng(3);
            n = 10; dur = 300;
            [cells, st, su] = test_AnalysisChoices.poissonTrains(n, dur, 2);
            % unit 2 fires 5 ms after half of unit 1's spikes
            extra = cells{1}(1:2:end) + 0.005;
            cells{2} = sort([cells{2}; extra]);
            st = [st; extra]; su = [su; 2 * ones(numel(extra), 1)];
            [st, o] = sort(st); su = su(o);
            p = struct('MaxLag', 0.05, 'SurrogateMethod', 'jitter', 'JitterSD', [], ...
                'N_Surrogates', 40, 'Percentile', 95, 'EdgeTest', 'fdr', 'FDR', 0.05);
            c = computeConnectivitySTTC(cells, st, su, n, dur, 20000, p);
            tc.verifyTrue(c.bu(1, 2) && c.bu(2, 1));
            tc.verifyLessThan(c.q(1, 2), 1e-3);
        end

        function testFeatureVersionStampedAndSaved(tc)
            proc = test_DataLayerFixes.processor(tc.Root, "ISIN");
            v = RecordingProcessor.FeatureDefinitionVersion;
            tc.verifyEqual(proc.FeatureVersion.UnitFeatures, v);
            f = fullfile(tc.Root, "p.mat");
            proc.save(f);
            loaded = RecordingProcessor.load(f);
            tc.verifyEqual(loaded.FeatureVersion.UnitFeatures, v);
        end

        function testMixedFeatureVersionsWarn(tc)
            [p1, p2] = test_DataLayerFixes.twoProcessors(tc.Root);
            warning('on', 'FeatureStore:mixedFeatureVersions');
            tc.verifyWarningFree(@() FeatureStore.checkFeatureVersions([p1 p2]));
            p2.FeatureVersion = struct();           % as saved before versions were stamped
            tc.verifyWarning(@() FeatureStore.checkFeatureVersions([p1 p2]), ...
                'FeatureStore:mixedFeatureVersions');
        end

        function testStratifiedEnrichmentIgnoresCultureRate(tc)
            % community 1 holds most units of culture A (responder rate 25%), community 2 most
            % of culture B (5%); within each culture responders are spread evenly over both
            % communities, so neither community is enriched once cultures are accounted for
            cult = [repmat("A", 200, 1); repmat("B", 200, 1)];
            M = [ones(160, 1); 2 * ones(40, 1); ones(40, 1); 2 * ones(160, 1)];
            resp = false(400, 1);
            resp([1:4:160, 161:4:200]) = true;              % 25% in each community of A
            resp([201:20:240, 241:20:400]) = true;          % 5% in each community of B
            elig = true(400, 1);
            [~, inh_pooled] = CellTypeClassifier.communityEnrichment(M, elig, resp, 0.05);
            [T, inh_strat]  = CellTypeClassifier.communityEnrichment(M, elig, resp, 0.05, cult);
            tc.verifyEqual(inh_pooled, 1);                  % pooled: confounded by culture
            tc.verifyEmpty(inh_strat);
            tc.verifyEqual(T.Expected, [42; 18], 'AbsTol', 1e-9);   % = observed counts
        end

        function testSingleStratumMatchesHypergeometric(tc)
            rng(4);
            M = randi(5, 300, 1); resp = rand(300, 1) < 0.2; elig = true(300, 1);
            T = CellTypeClassifier.communityEnrichment(M, elig, resp, 0.05, ones(300, 1));
            n_c = accumarray(M, 1); r_c = accumarray(M(resp), 1, [5 1]);
            p_ref = hygecdf(r_c - 1, 300, sum(resp), n_c, 'upper');
            tc.verifyEqual(T.pEnriched, p_ref, 'AbsTol', 1e-10);
            p_dep = hygecdf(r_c, 300, sum(resp), n_c);
            tc.verifyEqual(T.pDepleted, p_dep, 'AbsTol', 1e-10);
        end

        function testLvRIndependentOfQCRefractory(tc)
            sp = make_mock_sorting(tc.Root, 6, 300, 20000, 120);
            sd = SpikeData.fromKilosort(sp, struct('ChipID', 'C1', 'PlatingDate', '2024-01-01'));
            lvr = cell(1, 2);
            refr = [0.001 0.003];
            for k = 1:2
                p = struct(); p.QC.N_Units = 1; p.QC.RPV = 1; p.QC.RefractoryPeriod = refr(k);
                p.Analyses.Catch22 = false;
                proc = RecordingProcessor(sd, p);
                proc.runQC(); proc.computeUnitFeatures();
                lvr{k} = proc.UnitFeatureTable.RevisedLocalVariation;
            end
            tc.verifyEqual(lvr{1}, lvr{2});
            d = RecordingProcessor.returnDefaultParams();
            tc.verifyEqual(d.Activity.LvRRefractory, 0.005);
        end

        function testBurstParticipationGroup(tc)
            tc.verifyEqual(FeatureCatalog.features("BurstParticipation", "full"), "BurstFractionSpikes");
            fs = FeatureStore();
            fs.UnitTable = table(["u1"; "u2"], ["r"; "r"], [0.2; 0.7], [1; 2], ...
                'VariableNames', {'UnitID', 'RecordingID', 'BurstFractionSpikes', 'FiringRate'});
            X = fs.unitMatrix("BurstParticipation");
            tc.verifyEqual(string(X.Properties.VariableNames), "BurstFractionSpikes");
        end


        function testSilentUnitHasPValueOne(tc)
            rng(5);
            pre = poissrnd(2, 10, 40); post = poissrnd(2, 10, 40);
            pre(3, :) = 0; post(3, :) = 0;                  % silent unit
            r = bootstrapFiringRateResponse(pre, post, 200, 1e-3, 'NetworkCorrection', true);
            tc.verifyEqual(r.p_values(3), 1);
            tc.verifyTrue(ismember(3, r.unchanged));
        end

        function testCellTypeBurstFeaturesAgainstNull(tc)
            [st, su, lab, b, dur] = test_burstEI.network(1);
            [N, U] = computeCellTypeBurstEI(st, su, lab, b, dur, NNull=50);
            [net_names, unit_names] = cellTypeBurstFeatureNames();
            tc.verifyEqual(string(N.Properties.VariableNames), net_names);
            tc.verifyEqual(string(U.Properties.VariableNames), unit_names);
            tc.verifyEqual(height(U), numel(lab));
            % I fires before onset and is silent in bursts: lower participation and earlier
            % latency than E, beyond what firing rate explains
            tc.verifyLessThan(N.BurstUnitParticipation_IminusE_Excess, -0.3);
            tc.verifyLessThan(N.BurstUnitLatency_IminusE_Excess, -0.03);
            tc.verifyLessThan(N.BurstEI_IEContrastIn_Excess, 0);
        end

        function testCellTypeBurstExcessNearZeroForRandomLabels(tc)
            [st, su, lab, b, dur] = test_burstEI.network(2);
            rs = RandStream('twister', 'Seed', 5);
            e = find(lab == 1); fake = nan(size(lab)); fake(e) = 1;
            fake(e(randperm(rs, numel(e), 20))) = 2;
            N = computeCellTypeBurstEI(st, su, fake, b, dur, NNull=100);
            tc.verifyLessThan(abs(N.BurstUnitParticipation_IminusE_Excess), 0.05);
            tc.verifyLessThan(abs(N.BurstUnitLatency_IminusE_Excess), 0.02);
        end

        function testCellTypeBurstTooFewBurstsIsNaN(tc)
            [st, su, lab, b, dur] = test_burstEI.network(3);
            b = struct('Onset', b.Onset(1:3), 'Offset', b.Offset(1:3), 'T_peak', b.T_peak(1:3));
            N = computeCellTypeBurstEI(st, su, lab, b, dur, NNull=10);
            tc.verifyTrue(all(isnan(N{1, :})));
        end

        function testProcessorWritesCellTypeBurstColumns(tc)
            proc = test_DataLayerFixes.processor(tc.Root, "ISIN");
            proc.computeNetworkFeatures();
            proc.computeConnectivity();
            labels = ones(1, numel(proc.Units)); labels(1:3) = 2;
            proc.CellTypeLabels = labels;
            proc.computeCellTypeFeatures();
            tc.verifyEqual(proc.Status.CellTypeFeatures, "done");
            [net_names, unit_names] = cellTypeBurstFeatureNames();
            tc.verifyTrue(all(ismember(net_names, string(proc.NetworkFeatureTable.Properties.VariableNames))));
            tc.verifyTrue(all(ismember(unit_names, string(proc.UnitFeatureTable.Properties.VariableNames))));
            tc.verifyFalse(any(startsWith(string(proc.NetworkFeatureTable.Properties.VariableNames), "BurstLead")));
        end

    end

    methods (Static)

        function [cells, st, su] = poissonTrains(n, dur, rate)
            cells = cell(1, n); st = []; su = [];
            for u = 1:n
                t = sort(rand(poissrnd(rate * dur), 1) * dur);
                cells{u} = t;
                st = [st; t]; su = [su; u * ones(numel(t), 1)]; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
        end

    end
end
