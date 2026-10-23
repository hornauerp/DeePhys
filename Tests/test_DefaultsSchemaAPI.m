classdef test_DefaultsSchemaAPI < matlab.unittest.TestCase
% TEST_DEFAULTSSCHEMAAPI  Processing defaults, the feature schema without cell-type labels or
% bursts, RegularityFit gating, FeatureCatalog group resolution, option checks, compact
% stored models, the CCG Fs option and small API errors.
%
% Run: runtests('test_DefaultsSchemaAPI')

    properties
        Procs   % two processed mock recordings (shared, read-only)
    end

    methods (TestClassSetup)
        function makeProcessors(tc)
            w = warning('off', 'all'); c = onCleanup(@() warning(w));
            root = string(tempname); mkdir(root);
            tc.addTeardown(@() rmdir(root, 's'));
            p = struct();
            p.QC.N_Units = 1; p.QC.Amplitude = []; p.QC.RPV = 1; p.Analyses.Catch22 = false;   % random mock trains violate the refractory period
            for k = 1:2
                ks = make_mock_sorting(fullfile(root, "r" + k), 8, 400, 20000, 60);
                t = readNPY(fullfile(ks, 'spike_times.npy'));
                writeNPY(t + uint64(k - 1), fullfile(ks, 'spike_times.npy'));   % distinct recordings (IDs)
                sd = SpikeData.fromKilosort(ks, struct('ChipID', "C" + k, 'PlatingDate', '240101', ...
                    'DIV', 14, 'Group', "G" + k));
                proc = RecordingProcessor(sd, p);
                evalc('proc.runAll();');
                procs(k) = proc; %#ok<AGROW>
            end
            tc.Procs = procs;
        end
    end

    methods (Test)

        function testDefaults(tc)
            p = RecordingProcessor.returnDefaultParams();
            tc.verifyEqual(string(p.Bursts.Method), "ISINJoint");
            tc.verifyEqual(p.QC.FiringRate, [0.01 100]);
            tc.verifyEqual(p.STTC.MaxLag, 0.01);
            tc.verifyFalse(isfield(p, 'Save'));
            tc.verifyEqual(tc.Procs(1).Status.ParentFeatures, "skipped");
        end

        function testLogISIDeprecated(tc)
            [st, su] = simulateNetworkBursts(Mode="none", Seed=1, Duration=60);
            n = max(su); fr = accumarray(su, 1, [n 1]) / 60; cv = ones(n, 1);
            tc.verifyWarning(@() computeBursts(st, su, fr, cv, n, 60, struct('Method', "logISI")), ...
                'computeBursts:logISIDeprecated');
        end

        function testNoLabelColumnsWithoutLabels(tc)
            for proc = tc.Procs
                names = [string(proc.NetworkFeatureTable.Properties.VariableNames), ...
                         string(proc.UnitFeatureTable.Properties.VariableNames)];
                tc.verifyFalse(any(endsWith(names, ["_E", "_I", "_EE", "_II", "_EI"])), ...
                    strjoin(names(endsWith(names, ["_E", "_I", "_EE", "_II", "_EI"])), ', '));
                tc.verifyFalse(any(ismember(names, ["FractionExcNN", "SpatialMixingIndex", "BurstOriginExcFraction"])));
                % spatial burst columns exist for every recording (NaN without bursts)
                tc.verifyTrue(ismember("MeanBurstPropagationSpeed", names));
            end
        end

        function testBurstFractionZeroWithoutBursts(tc)
            for proc = tc.Procs
                if isempty(proc.Bursts.T_start)
                    tc.verifyEqual(proc.UnitFeatureTable.BurstFractionSpikes, zeros(numel(proc.Units), 1));
                end
            end
        end

        function testRegularityFitNeedsRhythm(tc)
            rng(2);
            noise = double(rand(1, 6000) < 0.05);                          % aperiodic
            [~, ~, fit_noise] = computeRegularity(noise, 0.1, 20);
            tc.verifyTrue(isnan(fit_noise));
            t = (0:5999) * 0.1;
            rhythm = double(mod(t, 2) < 0.2) + 0.05 * rand(1, 6000);       % 0.5 Hz bursts
            [f, ~, fit_rhythm] = computeRegularity(rhythm, 0.1, 20);
            tc.verifyEqual(f, 0.5, 'AbsTol', 0.02);
            tc.verifyTrue(isfinite(fit_rhythm));
            tc.verifyFalse(ismember("RegularityFit", FeatureCatalog.features("all", "core")));
        end

        function testCatalogResolution(tc)
            cols = ["Density_CCG", "Density_STTC_w", "RichClub1_CCG", "Density_CCG_EE", ...
                    "ConvexHullArea", "ConvexHullArea_E", "SpatialDecayTau_EE", "Waveform12", ...
                    "WaveformX", "ACG3", "FiringRate"];
            tc.verifyEqual(FeatureCatalog.columns("NetworkGraph", cols), ["Density_CCG", "Density_STTC_w", "RichClub1_CCG"]);
            tc.verifyEqual(FeatureCatalog.columns("CellTypeGraph", cols), "Density_CCG_EE");
            tc.verifyEqual(FeatureCatalog.columns("SpatialFeatures", cols), ["ConvexHullArea", "ConvexHullArea_E", "SpatialDecayTau_EE"]);
            tc.verifyEqual(FeatureCatalog.columns("ReferenceWaveform", cols), "Waveform12");
            tc.verifyEqual(FeatureCatalog.columns("ACG", cols), "ACG3");
            % core = curated names; the E/I variants of a core spatial name are not core
            core = cols(FeatureCatalog.inSet(cols, "core"));
            tc.verifyTrue(ismember("ConvexHullArea", core));
            tc.verifyFalse(any(ismember(["ConvexHullArea_E", "SpatialDecayTau_EE"], core)));
            tc.verifyTrue(ismember("Density_CCG", core));
            % every core name lies in its group's full set
            for g = FeatureCatalog.groups()
                c = FeatureCatalog.features(g, "core");
                tc.verifyTrue(all(FeatureCatalog.matchColumns(g, c) | ismember(c, FeatureCatalog.features(g, "full"))), g);
            end
        end

        function testFeatureStoreGroups(tc)
            fs = FeatureStore.fromProcessors(tc.Procs);
            tc.verifyError(@() fs.recordingMatrix("NotAGroup"), 'FeatureStore:unknownFeatureGroup');
            Xd = fs.recordingMatrix("Density");                             % column prefix still works
            tc.verifyTrue(all(startsWith(string(Xd.Properties.VariableNames), "Density")));
            Xg = fs.recordingMatrix("NetworkGraph");
            tc.verifyNotEmpty(Xg.Properties.VariableNames);
            report = evalc('r = fs.summarizeFeatureGroups();'); %#ok<NASGU>
            r = fs.summarizeFeatureGroups();
            tc.verifyEqual(r.recording.NetworkGraph, width(Xg));
            tc.verifyEqual(r.unit.ActivityFeatures, width(fs.unitMatrix("ActivityFeatures")));
        end

        function testUnknownOptionWarns(tc)
            [X, Y] = test_DefaultsSchemaAPI.toyData();
            o = struct('Algorithm', "svm", 'KFold', 3, 'ComputeImportance', false);
            o.Kfold = 2;                                                    % typo
            tc.verifyWarning(@() Classifier.classify(X, Y, o), 'Classifier:unknownOption');
            tc.verifyWarning(@() Regressor.regress(X, double(Y == "B"), setfield(o, 'Algorithm', "svm")), ...
                'Regressor:unknownOption'); %#ok<SFLD>
        end

        function testRegressorRFOptionsAndCompactModels(tc)
            [X, Y] = test_DefaultsSchemaAPI.toyData();
            y = double(Y == "B") + 0.1 * randn(numel(Y), 1);
            w = warning('off', 'all'); c = onCleanup(@() warning(w));
            tc.verifyError(@() Regressor.regress(X, y, struct('RF', struct('MinLeaf', 1))), 'Regressor:rfOptions');
            r = Regressor.regress(X, y, struct('KFold', 3, 'Seed', 1, 'RF', struct('MinLeafSize', 1, 'MinParentSize', 2)));
            tc.verifyTrue(contains(class(r(1).Mdl), "Compact"));
            rc = Classifier.classify(X, Y, struct('KFold', 3, 'Seed', 1, 'ComputeImportance', false));
            tc.verifyTrue(contains(class(rc(1).Mdl), "Compact"));
            tc.verifyEqual(Classifier.predict(rc(1), X), Classifier.predict(rc(1), X));   % compact model predicts
            rf = Classifier.classify(X, Y, struct('KFold', 3, 'Seed', 1, 'ComputeImportance', false, 'KeepFullModel', true));
            tc.verifyFalse(contains(class(rf(1).Mdl), "Compact"));
        end

        function testCCGHonoursFs(tc)
            % the option was never parsed (case mismatch); an invalid value now errors
            tc.verifyError(@() CCG([0.1; 0.2], [1; 1], 'binSize', 0.001, 'duration', 0.02, 'Fs', -1), ?MException);
            [a, t] = CCG([0.1; 0.1 + 0.0005], [1; 1], 'binSize', 0.0005, 'duration', 0.004, 'Fs', 1 / 30000);
            tc.verifyEqual(a(t == 0.0005), 1);
        end

        function testCultureLevelWithoutCompleteCultures(tc)
            w = warning('off', 'all'); c = onCleanup(@() warning(w));
            exp = Experiment.fromProcessors(tc.Procs);
            tc.verifyError(@() exp.classify("Culture", "Group", struct()), 'Experiment:noCultures');
        end

    end

    methods (Static)
        function [X, Y] = toyData()
            rng(3);
            Y = repelem(["A"; "B"], 15);
            X = table(double(Y == "B") * 3 + randn(30, 1), randn(30, 1), 'VariableNames', ["f1", "f2"]);
        end
    end
end
