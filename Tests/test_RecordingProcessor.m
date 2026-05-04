classdef test_RecordingProcessor < matlab.unittest.TestCase
% TEST_RECORDINGPROCESSOR  Integration tests for RecordingProcessor pipeline.
%
% Uses make_mock_sorting to build a synthetic Kilosort output directory,
% then exercises the full processing chain without external toolbox deps.
%
% Run: runtests('test_RecordingProcessor')

    methods (Test)

        function testConstructionFromKilosort(tc)
            [proc, sorting_path] = test_RecordingProcessor.buildTestProcessor();
            cleanup = onCleanup(@() rmdir(sorting_path, 's'));
            tc.verifyNotEmpty(proc.SpikeData);
            tc.verifyNotEmpty(proc.SpikeData.RecordingID);
            tc.verifyTrue(isstruct(proc.Parameters));
            % All status fields should start as "pending"
            fields = fieldnames(proc.Status);
            for f = 1:numel(fields)
                tc.verifyEqual(proc.Status.(fields{f}), "pending");
            end
        end

        function testRunQCPopulatesUnits(tc)
            [proc, sorting_path] = test_RecordingProcessor.buildTestProcessor();
            cleanup = onCleanup(@() rmdir(sorting_path, 's'));
            proc.runQC();
            tc.verifyEqual(proc.Status.QC, "done");
            tc.verifyGreaterThanOrEqual(numel(proc.Units), 1);
            for i = 1:numel(proc.Units)
                tc.verifyFalse(isempty(proc.Units(i).SpikeTimes));
                tc.verifyFalse(isempty(proc.Units(i).ReferenceWaveform));
                tc.verifyFalse(isempty(proc.Units(i).UnitID));
            end
        end

        function testComputeUnitFeatures(tc)
            [proc, sorting_path] = test_RecordingProcessor.buildTestProcessor();
            cleanup = onCleanup(@() rmdir(sorting_path, 's'));
            proc.runQC();
            proc.computeUnitFeatures();
            tc.verifyEqual(proc.Status.UnitFeatures, "done");
            tc.verifyFalse(isempty(proc.UnitFeatureTable));
            tc.verifyEqual(height(proc.UnitFeatureTable), numel(proc.Units));
            expected_cols = ["FiringRate", "MeanInterSpikeInterval", "HalfWidth"];
            col_names = string(proc.UnitFeatureTable.Properties.VariableNames);
            for c = 1:numel(expected_cols)
                tc.verifyTrue(ismember(expected_cols(c), col_names), ...
                    "Expected column missing: " + expected_cols(c));
            end
        end

        function testComputeNetworkFeatures(tc)
            [proc, sorting_path] = test_RecordingProcessor.buildTestProcessor();
            cleanup = onCleanup(@() rmdir(sorting_path, 's'));
            proc.runQC();
            proc.computeUnitFeatures();
            proc.computeNetworkFeatures();
            tc.verifyEqual(proc.Status.NetworkFeatures, "done");
            tc.verifyTrue(istable(proc.NetworkFeatureTable));
        end

        function testComputeConnectivity(tc)
            [proc, sorting_path] = test_RecordingProcessor.buildTestProcessor();
            cleanup = onCleanup(@() rmdir(sorting_path, 's'));
            proc.runQC();
            proc.computeConnectivity();
            tc.verifyEqual(proc.Status.Connectivity, "done");
            if numel(proc.Units) >= 2
                tc.verifyTrue(isfield(proc.Connectivity, 'CCG'));
            end
        end

        function testFullPipelineSequence(tc)
            [proc, sorting_path] = test_RecordingProcessor.buildTestProcessor();
            cleanup = onCleanup(@() rmdir(sorting_path, 's'));
            proc.runQC();
            proc.computeUnitFeatures();
            proc.computeNetworkFeatures();
            proc.computeConnectivity();
            tc.verifyEqual(proc.Status.QC, "done");
            tc.verifyEqual(proc.Status.UnitFeatures, "done");
            tc.verifyEqual(proc.Status.NetworkFeatures, "done");
            tc.verifyEqual(proc.Status.Connectivity, "done");
        end

        function testSaveLoadRoundTrip(tc)
            [proc, sorting_path] = test_RecordingProcessor.buildTestProcessor();
            cleanup_sort = onCleanup(@() rmdir(sorting_path, 's'));
            proc.runQC();
            proc.computeUnitFeatures();

            tmp = [tempname '.mat'];
            cleanup_mat = onCleanup(@() delete(tmp));
            proc.save(tmp);
            proc2 = RecordingProcessor.load(tmp);

            tc.verifyEqual(proc2.Status.QC, "done");
            tc.verifyEqual(proc2.Status.UnitFeatures, "done");
            tc.verifyEqual(proc2.SpikeData.RecordingID, proc.SpikeData.RecordingID);
            tc.verifyEqual(numel(proc2.Units), numel(proc.Units));
            tc.verifyEqual(height(proc2.UnitFeatureTable), height(proc.UnitFeatureTable));
        end

        function testIdempotentRerun(tc)
            [proc, sorting_path] = test_RecordingProcessor.buildTestProcessor();
            cleanup = onCleanup(@() rmdir(sorting_path, 's'));
            proc.runQC();
            n_units_first = numel(proc.Units);
            proc.runQC();  % second call — must short-circuit
            tc.verifyEqual(numel(proc.Units), n_units_first);
            tc.verifyEqual(proc.Status.QC, "done");

            proc.computeUnitFeatures();
            h_first = height(proc.UnitFeatureTable);
            proc.computeUnitFeatures();  % second call — must short-circuit
            tc.verifyEqual(height(proc.UnitFeatureTable), h_first);
        end

    end

    methods (Static, Access = private)

        function [proc, sorting_path] = buildTestProcessor()
            sorting_path = make_mock_sorting(string(tempdir), 5, 50, 20000, 60);
            metadata = struct('ChipID', 'TestChip', 'PlatingDate', '2024-01-01');
            sd = SpikeData.fromKilosort(sorting_path, metadata);
            p.Analyses.Regularity   = false;
            p.Analyses.Catch22      = false;
            p.Analyses.Bursts       = true;
            p.Analyses.Connectivity = "CCG";
            p.GraphFeatures.SmallWorldnessIterations = 10;
            proc = RecordingProcessor(sd, p);
        end

    end
end
