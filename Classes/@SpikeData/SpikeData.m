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
        RecordingID         string  % Hash of the spike data (fromKilosort) or of the stored
                                    % InputPath (legacy objects); independent of mount point
        SortingID           string  % Identity of the spike sorting the templates belong to:
                                    % segments of one sorting share it, separately sorted
                                    % recordings (e.g. other days) differ. Part of the UnitID.
        InputPath           string  % Absolute path to Kilosort output directory
        ParentPath          string  % Path to parent (concatenated) Kilosort directory, "" if none
        Metadata            struct  % Recording metadata (ChipID, PlatingDate, DIV, Mutation, ...)
        SamplingRate        double  % Samples per second
        Duration            double  % Recording duration in seconds
        ElectrodeCoordinates double % (N_channels x 2) XY electrode positions
        SpikeTimes          double  % (N_spikes x 1) spike times in seconds
        SpikeUnits          double  % (N_spikes x 1) 1-indexed template assignment
        TemplateWaveforms   double  % (N_templates x N_samples x N_channels)
    end

    methods (Static)

        function sd = fromKilosort(input_path, metadata)
        % FROMKILOSORT  Load spike data from a Kilosort output directory.
        %
        %   sd = SpikeData.fromKilosort(input_path)
        %   sd = SpikeData.fromKilosort(input_path, metadata_struct)
        %
        % Reads spike_times.npy, spike_templates.npy (or spike_clusters.npy), templates.npy,
        % channel_positions.npy (optional) and params.py from input_path: Kilosort output or
        % any sorter exported to phy format (e.g. SpikeInterface export_to_phy). Bombcell /
        % Kilosort label tables are optional (RecordingProcessor QC.Bombcell, QC.KSLabel).
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

            % Electrode coordinates (only the spatial features need them)
            pos_file = fullfile(input_path, 'channel_positions.npy');
            if isfile(pos_file)
                sd.ElectrodeCoordinates = readNPY(pos_file);
            else
                warning('SpikeData:noChannelPositions', ...
                    '%s not found: spatial features will be unavailable.', pos_file);
            end

            % Spike times: raw sample indices → seconds
            raw_times     = readNPY(SpikeData.requireFile(input_path, 'spike_times.npy'));
            sd.SpikeTimes = double(raw_times(:)) / sd.SamplingRate;

            % Spike → template (row of templates.npy), 0-indexed in file → 1-indexed.
            % Kilosort writes spike_templates.npy; other sorters exported to phy format
            % (e.g. SpikeInterface export_to_phy) may only provide spike_clusters.npy,
            % whose cluster indices then address templates.npy as well.
            unit_file = fullfile(input_path, 'spike_templates.npy');
            if ~isfile(unit_file)
                unit_file = SpikeData.requireFile(input_path, 'spike_clusters.npy');
            end
            raw_units     = readNPY(unit_file);
            sd.SpikeUnits = double(raw_units(:)) + 1;
            assert(numel(sd.SpikeUnits) == numel(sd.SpikeTimes), 'SpikeData:spikeCount', ...
                '%s has %d entries but spike_times.npy has %d.', unit_file, ...
                numel(sd.SpikeUnits), numel(sd.SpikeTimes));

            % Template waveforms (N_templates x N_samples x N_channels)
            sd.TemplateWaveforms = readNPY(SpikeData.requireFile(input_path, 'templates.npy'));
            assert(isempty(sd.SpikeUnits) || max(sd.SpikeUnits) <= size(sd.TemplateWaveforms, 1), ...
                'SpikeData:templateIndex', ['%s refers to template %d, but templates.npy has %d. ' ...
                'Export one template per unit (e.g. SpikeInterface export_to_phy).'], ...
                unit_file, max(sd.SpikeUnits) - 1, size(sd.TemplateWaveforms, 1));

            % Duration from params.py n_samples (more accurate than max spike time,
            % which underestimates when the recording goes silent at the end).
            n_samples = SpikeData.parseNSamples(input_path);
            if ~isnan(n_samples)
                sd.Duration = n_samples / sd.SamplingRate;
            else
                sd.Duration = ceil(max(sd.SpikeTimes));
            end

            % Recording ID from the spike data, so the same sorting gets the same ID on
            % every machine and mount point (the path differs between Windows and Linux mounts)
            sd.RecordingID = SpikeData.contentID(raw_times, raw_units, sd.SamplingRate);
            sd.SortingID   = SpikeData.resolveSortingID(input_path, sd.ParentPath, sd.RecordingID);
        end

        function id = resolveSortingID(input_path, parent_path, own_id)
        % RESOLVESORTINGID  Sorting identity of a recording (see the SortingID property).
        %   A segment of a longer sorting (explicit/derived ParentPath, or a "segment_N" or
        %   split_sortings window folder inside a sorter output) gets the fingerprint of that
        %   parent sorting, so
        %   all its segments share it; any other recording is its own sorting (own_id).
            src = "";
            if strlength(string(parent_path)) > 0 && isfolder(parent_path)
                src = string(parent_path);
            else
                [up, leaf] = fileparts(char(input_path));
                is_window = ~isempty(regexp(leaf, '^segment_\d+$', 'once')) || ...
                    isfile(fullfile(char(input_path), 'split_window.txt'));   % written by split_sortings
                if is_window && isfile(fullfile(up, 'spike_times.npy'))
                    src = string(up);
                end
            end
            if src == ""
                id = string(own_id);
            else
                id = SpikeData.sortingFingerprint(src);
            end
        end

        function id = sortingFingerprint(folder)
        % SORTINGFINGERPRINT  Path-independent hash of a sorter output folder: size, first and
        %   last 4 kB of spike_times.npy and of spike_templates.npy (or spike_clusters.npy).
            fp = struct();
            names = ["spike_times.npy", "spike_templates.npy", "spike_clusters.npy"];
            for k = 1:numel(names)
                f = fullfile(folder, names(k));
                if ~isfile(f), continue; end
                info = dir(f);
                fid = fopen(f, 'r'); c = onCleanup(@() fclose(fid));
                head = fread(fid, min(4096, info.bytes), '*uint8');
                fseek(fid, -min(4096, info.bytes), 'eof');
                tail = fread(fid, min(4096, info.bytes), '*uint8');
                clear c
                fp.(matlab.lang.makeValidName(names(k))) = struct('Bytes', info.bytes, 'Head', head', 'Tail', tail');
            end
            assert(~isempty(fieldnames(fp)), 'SpikeData:noSorting', 'No spike_times.npy in %s', folder);
            id = string(paramHash(fp));
        end

        function key = legacySortingKey(input_path)
        % LEGACYSORTINGKEY  Sorting identity for legacy objects, which keep no spike files:
        %   the sorting folder (the parent of a "segment_N" folder, otherwise the input
        %   folder) as its last six path components, the same on every mount point.
            p = strrep(char(input_path), '\', '/');
            p = regexprep(p, '/+$', '');
            parts = strsplit(p, '/');
            if ~isempty(parts) && ~isempty(regexp(parts{end}, '^segment_\d+$', 'once'))
                parts(end) = [];
            end
            parts = parts(~cellfun(@isempty, parts));
            key = string(strjoin(parts(max(1, end - 5):end), '/'));
        end

        function id = contentID(raw_times, raw_units, sampling_rate)
        % CONTENTID  Path-independent recording ID: hash of a fingerprint of the spike data
        % (spike count, exact sum of the sample indices, spikes per template, and the
        % first and last 10,000 spikes). Identical sortings give identical IDs.
            t = uint64(raw_times(:));
            u = uint64(raw_units(:));
            k = min(10000, numel(t));
            counts = [];
            if ~isempty(u), counts = accumarray(double(u) + 1, 1)'; end
            fp = struct('SamplingRate', double(sampling_rate), 'NSpikes', numel(t), ...
                'SampleSum', sum(t, 'native'), 'Counts', counts, ...
                'Head', [t(1:k), u(1:k)], 'Tail', [t(end-k+1:end), u(end-k+1:end)]);
            id = string(paramHash(fp));
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
            sd.SortingID            = string(paramHash(struct('Sorting', char(SpikeData.legacySortingKey(sd.InputPath)))));
        end

    end

    methods (Static, Access = private)

        function sr = parseSamplingRate(input_path)
        % sample_rate from params.py (any spacing, CRLF or LF line ends)
            sr = SpikeData.paramsValue(input_path, 'sample_rate');
            if isnan(sr)
                error('SpikeData:noSampleRate', 'No sample_rate in %s.', ...
                    fullfile(input_path, 'params.py'));
            end
        end

        function n = parseNSamples(input_path)
        % n_samples (or n_samples_dat) from params.py; NaN when absent.
            n = SpikeData.paramsValue(input_path, 'n_samples\w*');
        end

        function v = paramsValue(input_path, key_pattern)
        % Numeric value of the first "key = value" line in params.py matching key_pattern.
            params_path = fullfile(input_path, 'params.py');
            if ~isfile(params_path)
                error('SpikeData:paramsNotFound', 'Cannot open %s', params_path);
            end
            tok = regexp(fileread(params_path), ['(?m)^\s*' key_pattern '\s*=\s*([-+0-9.eE]+)'], ...
                'tokens', 'once');
            v = NaN;
            if ~isempty(tok), v = str2double(tok{1}); end
        end

        function f = requireFile(input_path, name)
            f = fullfile(input_path, name);
            if ~isfile(f)
                error('SpikeData:missingFile', ['%s not found. Expected a Kilosort / phy-format ' ...
                    'folder: spike_times.npy, spike_templates.npy (or spike_clusters.npy), ' ...
                    'templates.npy, params.py.'], f);
            end
        end

    end
end
