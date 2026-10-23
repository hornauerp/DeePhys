classdef test_SorterFormats < matlab.unittest.TestCase
% TEST_SORTERFORMATS  Loading works for Kilosort and other phy-format sorters, with or
% without Bombcell / Kilosort labels, and the RecordingID does not depend on the path.
%
% Run: runtests('test_SorterFormats')

    properties
        Root
    end

    methods (TestMethodSetup)
        function makeRoot(tc)
            tc.Root = string(tempname);
            mkdir(tc.Root);
            tc.addTeardown(@() rmdir(tc.Root, 's'));
        end
    end

    methods (Test)

        function testPhyExportWithoutKilosortExtras(tc)
            % SpikeInterface-style export: spike_clusters.npy only, "key=value" params.py
            % with CRLF line ends and no n_samples, no channel_positions.npy.
            ks = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            ref = SpikeData.fromKilosort(ks);
            copyfile(fullfile(ks, 'spike_templates.npy'), fullfile(ks, 'spike_clusters.npy'));
            delete(fullfile(ks, 'spike_templates.npy'));
            delete(fullfile(ks, 'channel_positions.npy'));
            fid = fopen(fullfile(ks, 'params.py'), 'w');
            fprintf(fid, 'dat_path = ''recording.dat''\r\nn_channels_dat=4\r\nsample_rate=20000.0\r\n');
            fclose(fid);
            sd = tc.verifyWarning(@() SpikeData.fromKilosort(ks), 'SpikeData:noChannelPositions');
            tc.verifyEqual(sd.SamplingRate, 20000);
            tc.verifyEqual(sd.SpikeUnits, ref.SpikeUnits);
            tc.verifyEqual(sd.SpikeTimes, ref.SpikeTimes);
            proc = RecordingProcessor(sd, struct('QC', struct('N_Units', 1)));
            proc.runQC();
            tc.verifyNotEmpty(proc.Units);
        end

        function testLabelFiltersWithoutLabelFiles(tc)
            % Bombcell / Kilosort label QC enabled but never run: warn, keep the standard QC.
            ks = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            sd = SpikeData.fromKilosort(ks);
            base = RecordingProcessor(sd, struct('QC', struct('N_Units', 1)));
            base.runQC();
            p.QC.N_Units = 1;
            p.QC.Bombcell.Enable = true;
            p.QC.KSLabel.Enable  = true;
            proc = RecordingProcessor(sd, p);
            tc.verifyWarning(@() proc.runQC(), 'RecordingProcessor:bombcell');
            tc.verifyEqual(numel(proc.Units), numel(base.Units));
        end

        function testBombcellTypeNames(tc)
            % AcceptedTypes may be given as Bombcell names ("GOOD") or as codes (1).
            ks = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            fid = fopen(fullfile(ks, 'cluster_bc_unitType.tsv'), 'w');
            fprintf(fid, 'cluster_id\tbc_unitType\n0\tGOOD\n1\tGOOD\n2\tGOOD\n3\tNOISE\n4\tMUA\n');
            fclose(fid);
            sd = SpikeData.fromKilosort(ks);
            p.QC.N_Units = 1;
            p.QC.Amplitude = [];                 % mock templates have no uV scale
            p.QC.Bombcell.Enable = true;
            p.QC.Bombcell.AcceptedTypes = "GOOD";
            by_name = RecordingProcessor(sd, p);
            tc.verifyWarningFree(@() by_name.runQC());
            p.QC.Bombcell.AcceptedTypes = 1;
            by_code = RecordingProcessor(sd, p);
            by_code.runQC();
            tc.verifyNotEmpty(by_name.Units);
            tc.verifyTrue(all([by_name.Units.TemplateID] <= 3));
            tc.verifyEqual([by_name.Units.TemplateID], [by_code.Units.TemplateID]);
        end

        function testAmplitudeCriterionSkippedWithoutScale(tc)
            % Without whitening_mat_inv.npy / amplitudes.npy the templates have no uV scale:
            % the amplitude criterion is skipped with a warning instead of rejecting units.
            ks = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            base = RecordingProcessor(SpikeData.fromKilosort(ks), struct('QC', struct('N_Units', 1, 'Amplitude', [])));
            base.runQC();
            w = readNPY(fullfile(ks, 'templates.npy'));
            writeNPY(single(w * 1e-4), fullfile(ks, 'templates.npy'));
            sd = SpikeData.fromKilosort(ks);
            tc.verifyEmpty(sd.TemplateAmplitudes);
            proc = RecordingProcessor(sd, struct('QC', struct('N_Units', 1)));
            tc.verifyWarning(@() proc.runQC(), 'RecordingProcessor:noAmplitudeScale');
            tc.verifyEqual(numel(proc.Units), numel(base.Units));
        end

        function testNaNTemplateRejected(tc)
            % A template without a waveform (NaN) is rejected; the other units pass.
            ks = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            base = RecordingProcessor(SpikeData.fromKilosort(ks), struct('QC', struct('N_Units', 1, 'Amplitude', [])));
            base.runQC();
            w = readNPY(fullfile(ks, 'templates.npy'));
            w(5, :, :) = NaN;
            writeNPY(single(w), fullfile(ks, 'templates.npy'));
            proc = RecordingProcessor(SpikeData.fromKilosort(ks), struct('QC', struct('N_Units', 1)));
            w_state = warning('off', 'RecordingProcessor:noAmplitudeScale');
            c = onCleanup(@() warning(w_state));
            proc.runQC();
            kept = setdiff([base.Units.TemplateID], 5);
            tc.verifyEqual(sort([proc.Units.TemplateID]), sort(kept));
        end

        function testKilosortAmplitudeInMicrovolts(tc)
            % trough(template x whitening_mat_inv) x median scaling / scaleproc x gain, used
            % by QC.Amplitude and written as the Amplitude unit feature
            ks = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            tmpl = double(readNPY(fullfile(ks, 'templates.npy')));
            su = double(readNPY(fullfile(ks, 'spike_templates.npy'))) + 1;
            rng(5);
            winv = 200 * eye(4) + 20 * (rand(4) - 0.5);
            scale = 10 + 5 * rand(numel(su), 1);
            writeNPY(winv, fullfile(ks, 'whitening_mat_inv.npy'));
            writeNPY(scale, fullfile(ks, 'amplitudes.npy'));
            fid = fopen(fullfile(fileparts(ks), 'spikeinterface_recording.json'), 'w');
            fprintf(fid, '{"properties": {"gain_to_uV": [6.3, 6.3, 6.3, 6.3]}}');
            fclose(fid);
            c = onCleanup(@() delete(fullfile(fileparts(ks), 'spikeinterface_recording.json')));
            expected = zeros(5, 1);
            for u = 1:5
                U = squeeze(tmpl(u, :, :)) * winv;
                expected(u) = -min(U(:)) * median(scale(su == u)) / 200;
            end
            sd = SpikeData.fromKilosort(ks);
            tc.verifyEqual(sd.TemplateAmplitudes, expected, 'RelTol', 1e-5);
            tc.verifyEqual(sd.GainToUV, 6.3);
            uV = expected * 6.3;
            [~, smallest] = min(uV);
            lim = [min(uV(uV > min(uV))) - 1, 1e5];          % excludes only the smallest unit
            proc = RecordingProcessor(sd, struct('QC', struct('N_Units', 1, 'Amplitude', lim)));
            proc.runQC();
            tc.verifyFalse(ismember(smallest, [proc.Units.TemplateID]));
            proc.computeUnitFeatures("WaveformFeatures");
            tc.verifyEqual(proc.UnitFeatureTable.Amplitude, uV([proc.Units.TemplateID]), 'RelTol', 1e-5);
            % Kilosort 4: mean waveforms in whitened space, amplitudes.npy is a PC-feature norm
            % (not a scaling) and there is no scaleproc: trough(template x whitening_mat_inv).
            % Its params contain Infinity, which jsondecode rejects.
            fid = fopen(fullfile(ks, 'spikeinterface_params.json'), 'w');
            fprintf(fid, '{"sorter_name": "kilosort4", "sorter_params": {"artifact_threshold": Infinity}}');
            fclose(fid);
            ks4_expected = arrayfun(@(u) -min(squeeze(tmpl(u, :, :)) * winv, [], 'all'), (1:5)');
            tc.verifyEqual(SpikeData.fromKilosort(ks).TemplateAmplitudes, ks4_expected, 'RelTol', 1e-5);
            % native Kilosort 4 output (no SpikeInterface metadata): ops.npy marks it
            delete(fullfile(ks, 'spikeinterface_params.json'));
            writeNPY(0, fullfile(ks, 'ops.npy'));
            tc.verifyEqual(SpikeData.fromKilosort(ks).TemplateAmplitudes, ks4_expected, 'RelTol', 1e-5);
            % another sorter (unknown scaling convention): no amplitude
            delete(fullfile(ks, 'ops.npy'));
            fid = fopen(fullfile(ks, 'spikeinterface_params.json'), 'w');
            fprintf(fid, '{"sorter_name": "mountainsort5", "sorter_params": {}}');
            fclose(fid);
            tc.verifyEmpty(SpikeData.fromKilosort(ks).TemplateAmplitudes);
        end

        function testRecordingIDIsPathIndependent(tc)
            ks = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            copy = fullfile(tc.Root, "elsewhere", "sorting");
            mkdir(fileparts(copy));
            copyfile(ks, copy);
            id1 = SpikeData.fromKilosort(ks).RecordingID;
            id2 = SpikeData.fromKilosort(copy).RecordingID;
            tc.verifyEqual(id1, id2);
            t = readNPY(fullfile(copy, 'spike_times.npy'));
            t(end) = t(end) + 1;                                     % one spike moved by a sample
            writeNPY(t, fullfile(copy, 'spike_times.npy'));
            tc.verifyNotEqual(SpikeData.fromKilosort(copy).RecordingID, id1);
        end

        function testMissingRequiredFile(tc)
            ks = make_mock_sorting(tc.Root, 5, 50, 20000, 60);
            delete(fullfile(ks, 'templates.npy'));
            tc.verifyError(@() SpikeData.fromKilosort(ks), 'SpikeData:missingFile');
        end

    end
end
