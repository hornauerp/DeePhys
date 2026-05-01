classdef RecordingProcessor < handle
% RECORDINGPROCESSOR  Lazy computation engine for one recording.
%
% Replaces MEArecording's god-constructor pattern with explicit, on-demand
% computation steps. Each compute method checks Status and short-circuits
% if already done.
%
% Save/load format is plain MATLAB types (no handle objects) so no
% saveobj/loadobj hacks are needed.
%
% USAGE:
%   % From new Kilosort output
%   sd   = SpikeData.fromKilosort('/path/to/ks', metadata);
%   proc = RecordingProcessor(sd);
%   proc.runQC();
%   proc.computeUnitFeatures();
%   proc.computeNetworkFeatures();
%   proc.computeConnectivity();
%   proc.save('/path/to/ks/RecordingProcessor.mat');
%
%   % Load saved processor
%   proc = RecordingProcessor.load('/path/to/ks/RecordingProcessor.mat');
%
%   % Migrate from old MEArecording .mat
%   proc = RecordingProcessor.fromLegacyMat('/path/to/ks/MEArecording.mat');
%
%   % Batch load
%   procs = RecordingProcessor.loadMany({path1, path2, ...});
%
%   % Assemble FeatureStore for group analysis
%   fs = FeatureStore.fromProcessors(procs);

    properties
        SpikeData           SpikeData   % Raw recording data (value class)
        Units               UnitData    % (1 x N) QC-passing units (value class array)
        Parameters          struct      % Analysis parameters (merged with defaults)
        UnitFeatureTable    table       % (N_units x F) per-unit feature table
        NetworkFeatureTable table       % (1 x F) network-level feature table
        Connectivity        struct      % CCG / STTC / DDC connectivity results
        Bursts              struct      % Burst detection results
        Status              struct      % Tracks which analyses have been run
    end

    % =====================================================================
    % Construction
    % =====================================================================
    methods

        function proc = RecordingProcessor(spike_data, parameters)
        % RECORDINGPROCESSOR  Construct without running any computation.
            arguments
                spike_data  SpikeData = SpikeData()
                parameters  struct    = struct()
            end
            if isempty(spike_data.RecordingID)
                return  % empty construction for array pre-allocation
            end
            proc.SpikeData  = spike_data;
            proc.Parameters = parseStructParameters(RecordingProcessor.returnDefaultParams(), parameters);
            proc.Status     = RecordingProcessor.emptyStatus();
        end

    end

    % =====================================================================
    % Computation (lazy — each step checks Status before running)
    % =====================================================================
    methods

        function runQC(proc)
        % RUNQC  Filter templates by quality criteria, populate proc.Units.
        %   After this step: proc.Units is a UnitData array of passing units.
            if proc.Status.QC == "done"
                return
            end
            sd = proc.SpikeData;
            p  = proc.Parameters;

            % Template waveform matrix + reference electrodes
            tmpl = sd.TemplateWaveforms;
            tmpl = tmpl(:, sum(tmpl, [1,3], 'omitnan') ~= 0, :);  % remove zero-padded channels
            [max_amps, max_idx] = min(tmpl, [], [2,3], 'linear');
            max_amps = abs(max_amps * p.QC.LSB);
            [~, ~, ref_elec] = ind2sub(size(tmpl), max_idx);

            % Peak waveform (normalised)
            n_t = size(tmpl, 1);
            n_s = size(tmpl, 2);
            s_idx  = repmat(1:n_s, n_t, 1);
            t_idx  = repmat((1:n_t)', 1, n_s);
            c_idx  = repmat(ref_elec(:), 1, n_s);
            lin_idx = sub2ind(size(tmpl), t_idx, s_idx, c_idx);
            peak_wf = tmpl(lin_idx)';                            % n_samples × n_templates
            mx = max(abs(peak_wf));
            mx(mx == 0) = 1;        % guard: all-zero template → keep zeros, not NaN/Inf
            norm_wf = peak_wf ./ mx;

            % Firing rates + spike times per unit
            st  = sd.SpikeTimes;
            su  = sd.SpikeUnits;
            N   = size(tmpl, 1);
            unit_spike_times = arrayfun(@(x) st(su == x), 1:N, 'un', 0);
            firing_rates     = cellfun(@numel, unit_spike_times) / sd.Duration;

            % QC checks
            if isempty(p.QC.GoodUnits)
                bad_power = RecordingProcessor.checkPower(norm_wf, p.QC.PowerCutoff);
                bad_amp   = RecordingProcessor.checkAmplitude(max_amps, p.QC.Amplitude);
                bad_fr    = RecordingProcessor.checkFiringRate(firing_rates, p.QC.FiringRate);
                bad_rpv   = RecordingProcessor.checkRPV(unit_spike_times, p.QC.RPV, p.QC.RefractoryPeriod);
                [is_axon, is_noisy] = RecordingProcessor.checkWaveform(norm_wf, p.QC.Axon, p.QC.Noise, p.QC.NoiseCutout, sd.SamplingRate);
                spike_counts = cellfun(@numel, unit_spike_times);
                min_n = p.QC.MinSpikeCount;
                bad_spike_count = spike_counts(:) < min_n;
                good = ~(bad_power | bad_amp | bad_fr | bad_rpv | is_axon | is_noisy | bad_spike_count);
            else
                good = p.QC.GoodUnits;
            end
            fprintf('RecordingProcessor QC: %d/%d units passed\n', sum(good), N);

            % Optional KS label filter
            use_ks = isfield(p.QC, 'KSLabel') && p.QC.KSLabel.Enable;
            if use_ks
                try
                    [ks_labels, ks_ids] = RecordingProcessor.loadKSLabels(sd.InputPath);
                    ks_pass = false(N, 1);
                    accepted = string(p.QC.KSLabel.AcceptedLabels);
                    for k = 1:numel(ks_ids)
                        tid = ks_ids(k) + 1;
                        if tid >= 1 && tid <= N
                            ks_pass(tid) = ismember(ks_labels(k), accepted);
                        end
                    end
                    good = good & ks_pass;
                catch ME
                    warning('RecordingProcessor:ksLabel', '%s', ME.message);
                end
            end

            % Build UnitData array for passing units
            good_idx = find(good);
            rec_id   = sd.RecordingID;
            ud(numel(good_idx)) = UnitData();
            for u = 1:numel(good_idx)
                tid = good_idx(u);
                ud(u) = UnitData(tid, ref_elec(tid), unit_spike_times{tid}, ...
                                 norm_wf(:, tid), rec_id, sd.SamplingRate, sd.Duration);
                ud(u).UnitID = RecordingProcessor.computeStableID(sd.Metadata, tid, ref_elec(tid));
                if use_ks
                    match = find(ks_ids == (tid - 1), 1);
                    if ~isempty(match)
                        ud(u).KSLabel = ks_labels(match);
                    end
                end
            end
            proc.Units  = ud;
            proc.Status.QC = "done";
        end

        function computeUnitFeatures(proc, feature_groups)
        % COMPUTEUNITFEATURES  Compute per-unit features via standalone functions.
            arguments
                proc
                feature_groups string = "all"
            end
            if proc.Status.UnitFeatures == "done"
                return
            end
            if proc.Status.QC ~= "done"
                proc.runQC();
            end
            if isempty(proc.Units)
                proc.UnitFeatureTable = table();
                proc.Status.UnitFeatures = "done";
                return
            end

            if isscalar(feature_groups) && feature_groups == "all"
                feat_groups = ["ActivityFeatures","WaveformFeatures","RegularityFeatures","Catch22"];
            else
                feat_groups = feature_groups;
            end

            n_units = numel(proc.Units);
            sr      = proc.SpikeData.SamplingRate;
            norm_wf = [proc.Units.ReferenceWaveform];
            feature_tables = {};

            % Waveform features
            if any(feat_groups == "WaveformFeatures")
                wf_feats = computeWaveformFeatures(norm_wf, sr);
                wf_tables = cellfun(@(s) struct2table(s), wf_feats, 'un', 0);
                feature_tables{end+1} = vertcat(wf_tables{:});
            end

            % Activity features (+ optional regularity + catch22)
            act_params.FanoBinWidth = 0.1;
            if isfield(proc.Parameters, 'Activity') && isfield(proc.Parameters.Activity, 'FanoBinWidth')
                act_params.FanoBinWidth = proc.Parameters.Activity.FanoBinWidth;
            end
            act_params.RefractoryPeriod = proc.Parameters.QC.RefractoryPeriod;
            if proc.Parameters.Analyses.Regularity
                act_params.Regularity = proc.Parameters.Regularity;
            end
            if proc.Parameters.Analyses.Catch22
                act_params.Catch22 = proc.Parameters.Catch22;
            end

            act_tables = cell(n_units, 1);
            reg_tables = cell(n_units, 1);
            c22_tables = cell(n_units, 1);
            for u = 1:n_units
                [act_tables{u}, reg_tables{u}, c22_tables{u}] = ...
                    computeActivityFeatures(proc.Units(u).SpikeTimes, ...
                        proc.Units(u).RecordingDuration, act_params);
            end
            if any(feat_groups == "ActivityFeatures")
                feature_tables{end+1} = vertcat(act_tables{:});
            end
            if any(feat_groups == "RegularityFeatures") && proc.Parameters.Analyses.Regularity
                feature_tables{end+1} = vertcat(reg_tables{:});
            end
            if any(feat_groups == "Catch22") && proc.Parameters.Analyses.Catch22
                feature_tables{end+1} = vertcat(c22_tables{:});
            end

            % Aligned waveforms
            aligned_wf = alignWaveformsFromData(norm_wf, sr);
            var_names = "Waveform" + (1:size(aligned_wf, 1));
            feature_tables{end+1} = array2table(aligned_wf', "VariableNames", var_names);

            proc.UnitFeatureTable = [feature_tables{:}];
            proc.Status.UnitFeatures = "done";
        end

        function computeParentFeatures(proc, acg_params)
        % COMPUTEPARENTFEATURES  Compute unit features from the parent (concatenated) recording.
        %
        % Currently computes Parent_ACG1..N: autocorrelograms from the full spike
        % train of the parent Kilosort sorting (more spikes = better ACG estimate).
        % Future parent-derived features (Parent_FiringRate, etc.) extend this method.
        %
        % Parent path is read from proc.SpikeData.ParentPath (set during SpikeData
        % construction from metadata.ParentInputPath or auto-derived). If no valid
        % parent path is available, a warning is printed and the method returns without
        % adding columns — child features serve as automatic fallback in FeatureStore.
        %
        % Results are stored as Parent_ACG1..N columns in proc.UnitFeatureTable.
        % Multiple sibling processors sharing the same parent share the CCG computation
        % via ParentSpikeLoader's persistent session cache.
            arguments
                proc
                acg_params struct = struct()
            end
            if proc.Status.ParentFeatures == "done"
                return
            end
            if proc.Status.QC ~= "done"
                proc.runQC();
            end

            % Resolve ACG parameters (match CellTypeClassifier defaults)
            defaults.BinSize = 0.0005;
            defaults.Lag     = 0.1;
            p = parseStructParameters(defaults, acg_params);

            parent_path = proc.SpikeData.ParentPath;
            if parent_path == "" || ~isfolder(parent_path)
                warning('RecordingProcessor:noParent', ...
                    ['No valid parent path for recording %s — ' ...
                     'parent features unavailable, child features used as fallback.'], ...
                    proc.SpikeData.InputPath);
                proc.Status.ParentFeatures = "done";
                return
            end

            if isempty(proc.Units)
                proc.Status.ParentFeatures = "done";
                return
            end

            % Load full-recording ACGs via shared cache (computed once per parent+params)
            try
                acg_map = ParentSpikeLoader.loadFullACGs( ...
                    parent_path, proc.SpikeData.SamplingRate, p);
            catch ME
                warning('RecordingProcessor:parentACGFailed', ...
                    'Parent ACG computation failed for %s: %s', parent_path, ME.message);
                proc.Status.ParentFeatures = "failed";
                return
            end

            % Extract ACG per unit by TemplateID lookup
            n_units = numel(proc.Units);
            % Determine n_bins from any entry in the map (all should be same)
            map_keys = keys(acg_map);
            if isempty(map_keys)
                proc.Status.ParentFeatures = "done";
                return
            end
            n_bins = numel(acg_map(map_keys{1}));

            full_acgs = NaN(n_bins, n_units);
            for u = 1:n_units
                tid_key = num2str(proc.Units(u).TemplateID);
                if acg_map.isKey(tid_key)
                    full_acgs(:, u) = acg_map(tid_key);
                end
                % If not found: column stays NaN (unit QC-rejected in parent or missing)
            end

            % Build table with Parent_ prefix columns (storage-only convention)
            acg_names   = "Parent_ACG" + (1:n_bins);
            parent_tbl  = array2table(full_acgs', 'VariableNames', acg_names);

            % ── Future parent features ────────────────────────────────────────────
            % Add more Parent_* columns here as needed, e.g.:
            %   parent_tbl.Parent_FiringRate = ...
            %   parent_tbl.Parent_MeanISI    = ...

            % Merge into UnitFeatureTable (remove any stale Parent_* columns first)
            if ~isempty(proc.UnitFeatureTable)
                existing = string(proc.UnitFeatureTable.Properties.VariableNames);
                stale    = startsWith(existing, "Parent_");
                if any(stale)
                    proc.UnitFeatureTable(:, stale) = [];
                end
                proc.UnitFeatureTable = [proc.UnitFeatureTable, parent_tbl];
            else
                proc.UnitFeatureTable = parent_tbl;
            end
            proc.Status.ParentFeatures = "done";
        end

        function computeNetworkFeatures(proc)
        % COMPUTENETWORKFEATURES  Compute burst stats, regularity, Catch22 at network level.
            if proc.Status.NetworkFeatures == "done"
                return
            end

            spike_times = proc.SpikeData.SpikeTimes;
            duration    = proc.SpikeData.Duration;
            net_struct  = struct();

            if proc.Parameters.Analyses.Regularity
                try
                    binned   = histcounts(spike_times, 'BinLimits', [0 duration], ...
                        'BinWidth', proc.Parameters.Regularity.Binning);
                    norm_act = binned / max(max(binned), 1);
                    [freq, mag, fit_coeff] = computeRegularity(norm_act, ...
                        proc.Parameters.Regularity.Binning, proc.Parameters.Regularity.N_peaks);
                    net_struct.Regularity = struct2table(struct( ...
                        'NetworkRegularityFrequency', freq, ...
                        'NetworkRegularityMagnitude', mag, ...
                        'NetworkRegularityFit', fit_coeff));
                catch ME
                    warning('RecordingProcessor:networkFeatures', 'Regularity failed: %s', ME.message);
                end
            end

            if proc.Parameters.Analyses.Catch22
                try
                    net_struct.Catch22 = computeCatch22(spike_times, duration, proc.Parameters.Catch22.BinSize);
                catch ME
                    warning('RecordingProcessor:networkFeatures', 'Catch22 failed: %s', ME.message);
                end
            end

            if proc.Parameters.Analyses.Bursts
                try
                    if proc.Status.UnitFeatures ~= "done"
                        proc.computeUnitFeatures();
                    end
                    burst_params.Method = "ISIN";
                    if isfield(proc.Parameters, 'Bursts') && isfield(proc.Parameters.Bursts, 'Method')
                        burst_params.Method = proc.Parameters.Bursts.Method;
                    end
                    if isfield(proc.Parameters, 'Bursts') && isfield(proc.Parameters.Bursts, 'MergeFactor')
                        burst_params.MergeFactor = proc.Parameters.Bursts.MergeFactor;
                    end
                    burst_params.Outlier = proc.Parameters.Outlier;
                    burst_params.Binning = 0.01;

                    fr  = proc.UnitFeatureTable.FiringRate;
                    cv  = proc.UnitFeatureTable.CVInterSpikeInterval;
                    [bursts, burst_features] = computeBursts(spike_times, ...
                        proc.SpikeData.SpikeUnits, fr, cv, ...
                        numel(proc.Units), duration, burst_params);
                    proc.Bursts = bursts;
                    net_struct.BurstFeatures = burst_features;

                    if ~isempty(bursts.T_start)
                        n_units = numel(proc.Units);
                        burst_frac = NaN(n_units, 1);
                        in_burst = false(size(spike_times));
                        for b = 1:length(bursts.T_start)
                            in_burst = in_burst | (spike_times >= bursts.T_start(b) & spike_times <= bursts.T_end(b));
                        end
                        for u = 1:n_units
                            mask = proc.SpikeData.SpikeUnits == proc.Units(u).TemplateID;
                            n_spikes = sum(mask);
                            if n_spikes > 0
                                burst_frac(u) = sum(mask & in_burst) / n_spikes;
                            end
                        end
                        proc.UnitFeatureTable.BurstFractionSpikes = burst_frac;
                    end
                catch ME
                    warning('RecordingProcessor:networkFeatures', 'Bursts failed: %s', ME.message);
                end
            end

            % Flatten struct to a one-row table
            if ~isempty(fieldnames(net_struct))
                proc.NetworkFeatureTable = RecordingProcessor.structToOneRowTable(net_struct);
            else
                proc.NetworkFeatureTable = table();
            end
            proc.Status.NetworkFeatures = "done";
        end

        function computeConnectivity(proc, methods)
        % COMPUTECONNECTIVITY  Compute CCG / STTC / DDC connectivity via standalone functions.
            arguments
                proc
                methods string = proc.Parameters.Analyses.Connectivity
            end
            if proc.Status.Connectivity == "done"
                return
            end
            if isempty(methods)
                proc.Status.Connectivity = "done";
                return
            end
            if proc.Status.QC ~= "done"
                proc.runQC();
            end

            st  = proc.SpikeData.SpikeTimes;
            su  = proc.SpikeData.SpikeUnits;
            sr  = proc.SpikeData.SamplingRate;
            n   = numel(proc.Units);
            dur = proc.SpikeData.Duration;

            for a = 1:numel(methods)
                try
                    switch methods(a)
                        case "CCG"
                            proc.Connectivity.CCG = computeConnectivityCCG(st, su, n, sr, proc.Parameters.CCG);
                        case "FullCCG"
                            proc.Connectivity.FullCCG = computeConnectivityCCG(st, su, n, sr, ...
                                struct('BinSize', 0.0005, 'Duration', 0.1, ...
                                       'Conv_w', proc.Parameters.CCG.Conv_w, ...
                                       'Alpha', proc.Parameters.CCG.Alpha));
                        case "STTC"
                            unit_st = arrayfun(@(u) u.SpikeTimes, proc.Units, 'un', 0);
                            proc.Connectivity.STTC = computeConnectivitySTTC( ...
                                unit_st, st, su, n, dur, sr, proc.Parameters.STTC);
                        case "DDC"
                            proc.Connectivity.DDC = computeConnectivityDDC(st, su, n, dur, proc.Parameters.DDC);
                    end
                catch ME
                    warning('RecordingProcessor:connectivity', '%s: %s', methods(a), ME.message);
                end
            end

            % Graph features
            if ~isempty(proc.Connectivity) && ~isempty(fieldnames(proc.Connectivity))
                try
                    gp = struct();
                    if isfield(proc.Parameters, 'GraphFeatures') && isfield(proc.Parameters.GraphFeatures, 'SmallWorldnessIterations')
                        gp.SmallWorldnessIterations = proc.Parameters.GraphFeatures.SmallWorldnessIterations;
                    end
                    [nw_graph, unit_graph] = computeGraphFeatures(proc.Connectivity, gp);
                    if ~isempty(proc.NetworkFeatureTable)
                        proc.NetworkFeatureTable = [proc.NetworkFeatureTable, nw_graph];
                    else
                        proc.NetworkFeatureTable = nw_graph;
                    end
                    if ~isempty(unit_graph) && ~isempty(proc.UnitFeatureTable)
                        proc.UnitFeatureTable = [proc.UnitFeatureTable, unit_graph];
                    end
                catch ME
                    warning('RecordingProcessor:graphFeatures', '%s', ME.message);
                end
            end
            proc.Status.Connectivity = "done";
        end

        function runAll(proc)
        % RUNALL  Convenience: run all enabled analyses in sequence.
            proc.runQC();
            if ~isempty(proc.Units)
                proc.computeUnitFeatures();
                proc.computeParentFeatures();
                proc.computeNetworkFeatures();
                proc.computeConnectivity();
            end
        end

    end

    % =====================================================================
    % Persistence
    % =====================================================================
    methods

        function save(proc, file_path)
        % SAVE  Write all data as plain MATLAB types (no handle objects).
        %   Load with RecordingProcessor.load(path).
            arguments
                proc      RecordingProcessor
                file_path (1,1) string
            end
            % Convert value-class objects to plain structs for clean serialization
            SpikeDataStruct       = RecordingProcessor.valueToStruct(proc.SpikeData);   %#ok<PROP>
            UnitsStructArray      = arrayfun(@RecordingProcessor.valueToStruct, proc.Units); %#ok<PROP>
            UnitFeatureTable      = proc.UnitFeatureTable;    %#ok<PROP>
            NetworkFeatureTable   = proc.NetworkFeatureTable; %#ok<PROP>
            Connectivity          = proc.Connectivity;        %#ok<PROP>
            Bursts                = proc.Bursts;              %#ok<PROP>
            Parameters            = proc.Parameters;          %#ok<PROP>
            Status                = proc.Status;              %#ok<PROP>
            builtin('save', file_path, 'SpikeDataStruct', 'UnitsStructArray', ...
                'UnitFeatureTable', 'NetworkFeatureTable', ...
                'Connectivity', 'Bursts', 'Parameters', 'Status');
        end

    end

    methods (Static)

        function proc = load(file_path)
        % LOAD  Restore a RecordingProcessor saved with proc.save().
            arguments
                file_path (1,1) string
            end
            s    = builtin('load', file_path);
            proc = RecordingProcessor();

            if isfield(s, 'SpikeDataStruct')
                proc.SpikeData = SpikeData.fromStruct(s.SpikeDataStruct);
            else
                % Old MEArecording format — attempt transparent migration
                try
                    proc = RecordingProcessor.fromLegacyMat(file_path);
                    fprintf('RecordingProcessor.load: auto-migrated legacy file %s\n', file_path);
                    return
                catch ME
                    warning('RecordingProcessor:load', ...
                        'No SpikeDataStruct in %s and legacy migration failed: %s', ...
                        file_path, ME.message);
                end
            end
            if isfield(s, 'UnitsStructArray') && ~isempty(s.UnitsStructArray)
                proc.Units = arrayfun(@UnitData.fromStruct, s.UnitsStructArray);
            end
            if isfield(s, 'UnitFeatureTable')
                proc.UnitFeatureTable = s.UnitFeatureTable;
            end
            if isfield(s, 'NetworkFeatureTable')
                proc.NetworkFeatureTable = s.NetworkFeatureTable;
            end
            if isfield(s, 'Connectivity')
                proc.Connectivity = s.Connectivity;
            end
            if isfield(s, 'Bursts')
                proc.Bursts = s.Bursts;
            end
            if isfield(s, 'Parameters')
                proc.Parameters = s.Parameters;
            else
                proc.Parameters = RecordingProcessor.returnDefaultParams();
            end
            if isfield(s, 'Status')
                proc.Status = s.Status;
                % Back-fill ParentFeatures flag if loaded from older save (field didn't exist)
                if ~isfield(proc.Status, 'ParentFeatures')
                    % Back-fill missing field for files saved before tri-state Status.
                    if ~isempty(proc.UnitFeatureTable) && ...
                            any(startsWith(string(proc.UnitFeatureTable.Properties.VariableNames), "Parent_"))
                        proc.Status.ParentFeatures = "done";
                    else
                        proc.Status.ParentFeatures = "pending";
                    end
                end
                % Upgrade legacy boolean Status fields to tri-state strings.
                proc.Status = RecordingProcessor.upgradeLegacyStatus(proc.Status);
            else
                proc.Status = RecordingProcessor.emptyStatus();
            end
        end

        function procs = loadMany(file_paths)
        % LOADMANY  Load multiple processors in parallel.
            if ischar(file_paths) || isstring(file_paths)
                file_paths = cellstr(file_paths);
            end
            N = numel(file_paths);
            procs(N) = RecordingProcessor();
            parfor i = 1:N
                try
                    procs(i) = RecordingProcessor.load(string(file_paths{i}));
                catch ME
                    warning('RecordingProcessor:loadMany', 'Failed to load %s: %s', file_paths{i}, ME.message);
                end
            end
        end

        function proc = fromKilosort(input_path, metadata, parameters)
        % FROMKILOSORT  Create a RecordingProcessor from Kilosort output.
        %   Does NOT run any analysis — call proc.runAll() or individual steps.
            arguments
                input_path  (1,1) string
                metadata    (1,1) struct = struct()
                parameters  (1,1) struct = struct()
            end
            sd   = SpikeData.fromKilosort(input_path, metadata);
            proc = RecordingProcessor(sd, parameters);
        end

        function proc = fromLegacyMat(mat_path, parameters)
        % FROMLEGACYMAT  Migrate an old MEArecording.mat to the new format.
        %
        % Loads the old object, extracts SpikeData + UnitData + feature tables from
        % its already-computed properties. Sets Status = "migrated" so no recomputation.
            arguments
                mat_path   (1,1) string
                parameters (1,1) struct = struct()
            end
            data = builtin('load', mat_path);
            % Find the MEArecording object in the saved variables
            fns  = fieldnames(data);
            obj  = [];
            for i = 1:numel(fns)
                if isa(data.(fns{i}), 'MEArecording')
                    obj = data.(fns{i});
                    break
                end
            end
            if isempty(obj)
                error('RecordingProcessor:fromLegacyMat', ...
                    'No MEArecording object found in %s', mat_path);
            end

            % Build SpikeData from old object
            sd   = SpikeData.fromLegacyMEArecording(obj);
            proc = RecordingProcessor(sd, parameters);

            % Extract UnitData from old Unit objects
            if ~isempty(obj.Units)
                proc.Units = UnitData.fromLegacyUnitArray(obj.Units, sd.RecordingID);
                proc.Status.QC = "done";
            end

            % Extract unit feature table
            if ~isempty(obj.Units)
                try
                    proc.UnitFeatureTable = FeatureAssembly.unitFeatures(obj, "all");
                    proc.Status.UnitFeatures = "done";
                catch ME
                    warning('RecordingProcessor:fromLegacyMat', 'Unit features: %s', ME.message);
                end
            end

            % Extract recording-level / network features
            if ~isempty(obj.NetworkFeatures)
                try
                    proc.NetworkFeatureTable = RecordingProcessor.structToOneRowTable(obj.NetworkFeatures);
                    proc.Status.NetworkFeatures = "done";
                catch ME
                    warning('RecordingProcessor:fromLegacyMat', 'Network features: %s', ME.message);
                end
            end

            % Connectivity
            if ~isempty(obj.Connectivity) && ~isempty(fieldnames(obj.Connectivity))
                proc.Connectivity = obj.Connectivity;
                proc.Status.Connectivity = "done";
            end

            % Bursts
            if ~isempty(obj.Bursts)
                proc.Bursts = obj.Bursts;
            end

            fprintf('Migrated %s (%d units)\n', mat_path, numel(proc.Units));
        end

        function params = returnDefaultParams()
            params = MEArecording.returnDefaultParams();
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Static, Access = private)

        function s = emptyStatus()
            s.QC              = "pending";
            s.UnitFeatures    = "pending";
            s.ParentFeatures  = "pending";
            s.NetworkFeatures = "pending";
            s.Connectivity    = "pending";
        end

        function s = upgradeLegacyStatus(s)
        % Convert boolean Status fields saved before tri-state to string form.
            fields = fieldnames(s);
            for fi = 1:numel(fields)
                v = s.(fields{fi});
                if islogical(v) || isnumeric(v)
                    s.(fields{fi}) = "pending";
                    if v
                        s.(fields{fi}) = "done";
                    end
                end
            end
        end

        function uid = computeStableID(metadata, template_id, ref_electrode)
        % Compute deterministic StableID from recording metadata + unit identity.
            chip_id     = '';
            plating_date = '';
            if isfield(metadata, 'ChipID'),      chip_id     = char(string(metadata.ChipID));     end
            if isfield(metadata, 'PlatingDate'), plating_date = char(string(metadata.PlatingDate)); end
            uid = paramHash(struct('ChipID', chip_id, 'PlatingDate', plating_date, ...
                'TemplateID', template_id, 'ReferenceElectrode', ref_electrode));
        end

        function [labels, ids] = loadKSLabels(input_path)
            candidates = { fullfile(input_path, 'cluster_KSLabel.tsv'), ...
                           fullfile(fileparts(input_path), 'cluster_KSLabel.tsv') };
            path = '';
            for c = 1:numel(candidates)
                if exist(candidates{c}, 'file') == 2
                    path = candidates{c};
                    break
                end
            end
            if isempty(path)
                error('DeePhys:noKSLabels', 'cluster_KSLabel.tsv not found in %s', input_path);
            end
            T      = readtable(path, 'FileType', 'text', 'Delimiter', '\t');
            ids    = T.cluster_id;
            labels = string(T.KSLabel);
        end

        % ---- QC checks (ported from MEArecording) ----

        function bad = checkPower(norm_wf, cutoff)
            n     = size(norm_wf, 1);
            y     = fft(norm_wf);
            power = abs(y).^2 / n;
            bad   = max(power)' > cutoff;
            fprintf('QC power: %d bad\n', sum(bad));
        end

        function bad = checkAmplitude(amps, range)
            if isempty(range) || any(isnan(range))
                bad = false(size(amps));
            else
                bad = amps < range(1) | amps > range(2) | isnan(amps);
            end
            fprintf('QC amplitude: %d bad\n', sum(bad));
        end

        function bad = checkFiringRate(fr, range)
            if isempty(range) || any(isnan(range))
                bad = false(size(fr));
            else
                bad = fr < range(1) | fr > range(2);
            end
            fprintf('QC firing rate: %d bad\n', sum(bad));
        end

        function bad = checkRPV(spike_times_cell, rpv_thresh, ref_period)
            if isempty(rpv_thresh) || isnan(rpv_thresh)
                bad = false(numel(spike_times_cell), 1);
            else
                ratio = cellfun(@(x) sum(diff(x) < ref_period) / max(numel(x)-1, 1), spike_times_cell);
                bad   = ratio(:) > rpv_thresh;
            end
            fprintf('QC RPV: %d bad\n', sum(bad));
        end

        function [is_axon, is_noisy] = checkWaveform(norm_wf, axon_thresh, noise_thresh, noise_cutout, sampling_rate)
            if nargin < 5 || isempty(sampling_rate)
                sampling_rate = 10000;  % Maxwell default; pass sd.SamplingRate for other systems
            end
            SAMPLES_PER_MS = sampling_rate / 1000;
            if isempty(axon_thresh) || isnan(axon_thresh)
                is_axon = false(size(norm_wf, 2), 1);
            else
                amp_ratio = -(max(norm_wf) ./ min(norm_wf))';
                is_axon   = amp_ratio > axon_thresh;
            end
            fprintf('QC axonal: %d\n', sum(is_axon));

            if isempty(noise_thresh) || isnan(noise_thresh)
                is_noisy = false(size(norm_wf, 2), 1);
            else
                [~, peak_idx] = min(norm_wf);
                avg_peak = median(peak_idx(peak_idx > 1));
                if isnan(avg_peak)
                    is_noisy = false(size(norm_wf, 2), 1);
                else
                    cutout = max(1, round(avg_peak) + noise_cutout(1)*SAMPLES_PER_MS) : ...
                             min(size(norm_wf, 1), round(avg_peak) + noise_cutout(2)*SAMPLES_PER_MS);
                    if length(cutout) < 3
                        is_noisy = false(size(norm_wf, 2), 1);
                        fprintf('QC noisy: %d\n', sum(is_noisy));
                        return
                    end
                    bad_peak = peak_idx < cutout(1) | peak_idx > cutout(end);
                    n_changes = sum(abs(diff(diff(norm_wf(cutout, :)) > 0)))';
                    is_noisy  = n_changes > noise_thresh | bad_peak';
                end
            end
            fprintf('QC noisy: %d\n', sum(is_noisy));
        end

        function T = structToOneRowTable(s)
        % Flatten a (possibly nested) struct of scalar values into a one-row table.
            T = table();
            fns = fieldnames(s);
            for i = 1:numel(fns)
                v = s.(fns{i});
                if istable(v)
                    T = [T, v]; %#ok<AGROW>
                elseif isstruct(v)
                    T = [T, RecordingProcessor.structToOneRowTable(v)]; %#ok<AGROW>
                elseif isnumeric(v) && ~isempty(v)
                    if isscalar(v)
                        T.(fns{i}) = v;
                    else
                        var_names = string(fns{i}) + (1:numel(v));
                        T = [T, array2table(v(:)', 'VariableNames', var_names)]; %#ok<AGROW>
                    end
                end
            end
        end

        function s = valueToStruct(obj)
        % Convert a value-class object to a plain struct for serialization.
            s = struct();
            props = properties(obj);
            for i = 1:numel(props)
                s.(props{i}) = obj.(props{i});
            end
        end

    end
end
