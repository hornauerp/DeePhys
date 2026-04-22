classdef SpikeData
% SPIKEDATA  Pure data container for Kilosort spike-sorted output.
%
% Value class — no computation, no side effects. Trivially serializable
% because it holds only plain arrays and structs (no handle objects).
%
% USAGE:
%   sd = SpikeData.fromKilosort('/path/to/ks_output', metadata_struct);
%   sd = SpikeData.fromStruct(s);   % reconstruct from saved plain struct

    properties
        RecordingID         string  % Deterministic hash of InputPath
        InputPath           string  % Absolute path to Kilosort output directory
        ParentPath          string  % Path to parent (concatenated) Kilosort directory, "" if none
        Metadata            struct  % Recording metadata (ChipID, PlatingDate, DIV, Mutation, ...)
        SamplingRate        double  % Samples per second
        Duration            double  % Recording duration in seconds
        ElectrodeCoordinates double % (N_channels x 2) XY electrode positions
        SpikeTimes          double  % (N_spikes x 1) spike times in seconds
        SpikeUnits          double  % (N_spikes x 1) 1-indexed template assignment
        TemplateWaveforms   double  % (N_samples x N_channels x N_templates)
    end

    methods (Static)

        function sd = fromKilosort(input_path, metadata)
        % FROMKILOSORT  Load spike data from a Kilosort output directory.
        %
        %   sd = SpikeData.fromKilosort(input_path)
        %   sd = SpikeData.fromKilosort(input_path, metadata_struct)
        %
        % Reads spike_times.npy, spike_templates.npy, templates.npy,
        % channel_positions.npy, and params.py from input_path.
            arguments
                input_path  (1,1) string
                metadata    (1,1) struct = struct()
            end

            sd = SpikeData();
            sd.InputPath = input_path;

            % Resolve parent path (concatenated recording that this segment was split from).
            % Priority: explicit metadata.ParentInputPath > auto-derived sibling qc_output > ""
            if isfield(metadata, 'ParentInputPath') && ~isempty(metadata.ParentInputPath)
                sd.ParentPath = string(metadata.ParentInputPath);
            else
                candidate = fullfile(fileparts(char(input_path)), 'qc_output');
                if isfolder(candidate) && isfile(fullfile(candidate, 'spike_times.npy'))
                    sd.ParentPath = string(candidate);
                else
                    sd.ParentPath = "";
                end
            end

            % Merge InputPath into metadata
            metadata.InputPath = char(input_path);

            % Compute DIV from dates if not already present
            if ~isfield(metadata, 'DIV') || isempty(metadata.DIV)
                if isfield(metadata, 'RecordingDate') && ~isempty(metadata.RecordingDate) && ...
                   isfield(metadata, 'PlatingDate')   && ~isempty(metadata.PlatingDate)
                    try
                        rec_dt  = parseDateRobust(metadata.RecordingDate);
                        plat_dt = parseDateRobust(metadata.PlatingDate);
                        if ~isnat(rec_dt) && ~isnat(plat_dt)
                            metadata.DIV = days(rec_dt - plat_dt);
                        end
                    catch
                    end
                end
            end
            sd.Metadata = metadata;

            % Sampling rate from params.py
            sd.SamplingRate = SpikeData.parseSamplingRate(input_path);

            % Electrode coordinates
            sd.ElectrodeCoordinates = readNPY(fullfile(input_path, 'channel_positions.npy'));

            % Spike times: raw sample indices → seconds
            raw_times    = double(readNPY(fullfile(input_path, 'spike_times.npy')));
            sd.SpikeTimes = raw_times / sd.SamplingRate;

            % Template assignments: 0-indexed in file → 1-indexed
            sd.SpikeUnits = double(readNPY(fullfile(input_path, 'spike_templates.npy'))) + 1;

            % Template waveforms (N_samples x N_channels x N_templates)
            sd.TemplateWaveforms = readNPY(fullfile(input_path, 'templates.npy'));

            % Duration from params.py n_samples (more accurate than max spike time,
            % which underestimates when the recording goes silent at the end).
            n_samples = SpikeData.parseNSamples(input_path);
            if ~isnan(n_samples)
                sd.Duration = n_samples / sd.SamplingRate;
            else
                sd.Duration = ceil(max(sd.SpikeTimes));
            end

            % Deterministic recording ID
            sd.RecordingID = paramHash(struct('InputPath', char(input_path)));
        end

        function sd = fromStruct(s)
        % FROMSTRUCT  Reconstruct SpikeData from a plain struct (e.g. loaded from .mat).
            sd = SpikeData();
            fns = fieldnames(s);
            for i = 1:numel(fns)
                try
                    sd.(fns{i}) = s.(fns{i});
                catch
                end
            end
        end

        function sd = fromLegacyMEArecording(obj)
        % FROMLEGACYMEARECORDING  Extract SpikeData from an old MEArecording object.
        %   Used by RecordingProcessor.fromLegacyMat.
            arguments
                obj MEArecording
            end
            sd = SpikeData();
            sd.InputPath            = string(obj.Metadata.InputPath);
            sd.Metadata             = obj.Metadata;
            sd.SamplingRate         = obj.RecordingInfo.SamplingRate;
            sd.Duration             = obj.RecordingInfo.Duration;
            sd.ElectrodeCoordinates = obj.RecordingInfo.ElectrodeCoordinates;
            sd.SpikeTimes           = obj.Spikes.Times;
            sd.SpikeUnits           = obj.Spikes.Units;
            sd.TemplateWaveforms    = [];   % not preserved in old .mat (freed after construction)
            sd.RecordingID          = paramHash(struct('InputPath', char(sd.InputPath)));
        end

    end

    methods (Static, Access = private)

        function sr = parseSamplingRate(input_path)
            params_path = fullfile(input_path, 'params.py');
            fid  = fopen(params_path, 'r');
            txt  = fscanf(fid, '%c');
            fclose(fid);
            lines   = strsplit(txt, newline);
            sr_line = lines{startsWith(lines, 'sample')};
            parts   = strsplit(sr_line, ' = ');
            sr      = str2double(parts{2});
        end

        function n = parseNSamples(input_path)
        % Parse n_samples (or n_samples_dat) from params.py. Returns NaN when absent.
            params_path = fullfile(input_path, 'params.py');
            fid = fopen(params_path, 'r');
            txt = fscanf(fid, '%c');
            fclose(fid);
            lines = strsplit(txt, newline);
            idx   = find(startsWith(lines, 'n_samples'), 1);
            if isempty(idx)
                n = NaN;
                return
            end
            parts = strsplit(lines{idx}, ' = ');
            n     = str2double(parts{2});
        end

    end
end
