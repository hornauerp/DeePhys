classdef test_ExperimentCV < matlab.unittest.TestCase
% TEST_EXPERIMENTCV  Cross-validation groups of Experiment.prepareML.
%
% Recordings of one culture (DIVs, segments) share their neurons, so by default
% (CVLevel "culture") all rows of a culture must land in the same CV group, at
% recording and at unit level.
%
% Run: runtests('test_ExperimentCV')

    methods (Test)

        function testRecordingLevelGroupsByCulture(tc)
            exp = test_ExperimentCV.makeExperiment();
            [~, ~, g] = exp.prepareML("Recording", "Genotype", test_ExperimentCV.opts("culture"));
            tc.verifyEqual(numel(unique(g)), 2);
            tc.verifyEqual(g(1), g(2));          % r1, r2: culture A
            tc.verifyEqual(g(3), g(4));          % r3, r4: culture B
            tc.verifyNotEqual(g(1), g(3));
        end

        function testUnitLevelGroupsByCulture(tc)
            exp = test_ExperimentCV.makeExperiment();
            [~, ~, g] = exp.prepareML("Unit", "Genotype", test_ExperimentCV.opts("culture"));
            cul = extractBefore(exp.FeatureStore.UnitTable.RecordingID, 2);   % "a" / "b"
            tc.verifyEqual(numel(unique(g)), 2);
            tc.verifyEqual(numel(unique(g(cul == "a"))), 1);
        end

        function testRecordingLevelOption(tc)
            exp = test_ExperimentCV.makeExperiment();
            [~, ~, g] = exp.prepareML("Unit", "Genotype", test_ExperimentCV.opts("recording"));
            tc.verifyEqual(numel(unique(g)), 4);
        end

        function testMissingIdentityKeysFallsBackToRecording(tc)
            exp = test_ExperimentCV.makeExperiment();
            o = test_ExperimentCV.opts("culture"); o.IdentityKeys = "Well";
            [~, ~, g] = tc.verifyWarning(@() exp.prepareML("Recording", "Genotype", o), ...
                'Experiment:noCultureKeys');
            tc.verifyEqual(numel(unique(g)), 4);   % one group per recording
        end

    end

    methods (Static, Access = private)

        function o = opts(level)
            o = struct('FeatureGroups', "all", 'ParentFeatures', string.empty, ...
                'FeatureSet', "full", 'MinFiringRate', 0, 'CVLevel', level, ...
                'IdentityKeys', ["ChipID", "PlatingDate"], 'NormalizationVar', "");
        end

        function exp = makeExperiment()
            % Two cultures (A, B) with two recordings each; two units per recording.
            rec  = ["a1"; "a2"; "b1"; "b2"];
            chip = ["A"; "A"; "B"; "B"];
            meta = table(rec, chip, [1; 1; 2; 2], ["WT"; "WT"; "KO"; "KO"], ...
                'VariableNames', {'RecordingID', 'ChipID', 'PlatingDate', 'Genotype'});
            rec_tbl = [meta, table([1; 2; 3; 4], 'VariableNames', {'BurstRate'})];
            u_rec = repelem(rec, 2);
            unit_tbl = table(["u1"; "u2"; "u1"; "u2"; "u3"; "u4"; "u3"; "u4"], u_rec, ...
                repelem(meta.Genotype, 2), (1:8)', ...
                'VariableNames', {'UnitID', 'RecordingID', 'Genotype', 'FiringRate'});
            fs = FeatureStore();
            fs.UnitTable = unit_tbl; fs.RecordingTable = rec_tbl; fs.MetadataTable = meta;
            exp = Experiment();
            exp.FeatureStore = fs;
        end

    end
end
