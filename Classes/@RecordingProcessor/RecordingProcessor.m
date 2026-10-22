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
        Connectivity        struct = struct()   % CCG / STTC / DDC connectivity results (fields added per method)
        Bursts              struct      % Burst detection results
        CellTypeLabels      double      % (1 x N_units) 1=exc, 2=inh, NaN=unclassified
        Status              struct      % Tracks which analyses have been run
        FeatureVersion      struct = struct()   % FeatureDefinitionVersion each step ran under (missing = unknown, older)
    end

    properties (Constant)
        % Version of the feature definitions; raise it whenever a computed feature changes
        % meaning under the same column name, so FeatureStore.fromProcessors can refuse to mix
        % them silently.
        %   (unstamped) definitions before version 2
        %   2  AUC_* in amplitude x ms; HalfWidth and T2P searched after the trough within 2 ms;
        %      RegularityMagnitude = spectral peak / mean; SpatialMixingIndex relative to random
        %      labels; SpatialEIMoransI on unit labels; RichClub with fixed columns;
        %      BurstFractionSpikes; STTC graph edges by surrogate z-test with BH-FDR;
        %      RevisedLocalVariation with R = Activity.LvRRefractory (5 ms, was the 2 ms QC
        %      refractory period); spatial firing-rate Moran's I on log10 rate; cell-type
        %      burst features from burst onset / offset against a rate-matched null
        %      (computeCellTypeBurstEI: BurstEI_*, BurstUnit*), replacing BurstLeadTime_*,
        %      BurstLeadFraction, MeanBurstFraction_*, IntraBurstRate_*, BurstParticipation*
        %   3  catch22 without the catch24 extras (SC_DN_Mean / SC_DN_Spread_Std and the
        %      network DN_Mean / DN_Spread_Std are no longer written); SmallWorldness from
        %      path lengths of reachable pairs (was NaN or 0 with any unreachable pair); QC of a
        %      window of a longer sorting on the parent's spike trains (QC.RateReference)
        FeatureDefinitionVersion = 3
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
            % (templates x samples x channels). Zero-padded channels need no removal: they never
            % hold the trough, and dropping columns would break the electrode index of ref_elec.
            tmpl = sd.TemplateWaveforms;
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
            unit_spike_times = RecordingProcessor.splitByUnit(st, su, N);
            firing_rates     = cellfun(@numel, unit_spike_times) / sd.Duration;

            % QC checks
            if isempty(p.QC.GoodUnits)
                bad_power = RecordingProcessor.checkPower(norm_wf, p.QC.PowerCutoff);
                bad_amp   = RecordingProcessor.checkAmplitude(max_amps, p.QC.Amplitude);
                if N > 0 && all(bad_amp) && all(max_amps(:) < p.QC.Amplitude(1))
                    % Kilosort writes whitened templates; without whitening_mat_inv.npy
                    % their scale is not in ADC units, so a microvolt range rejects everything.
                    warning('RecordingProcessor:amplitudeScale', ['%s: all %d template amplitudes are below ' ...
                        'QC.Amplitude(1) = %g (max %.3g); the templates are probably whitened, so the ' ...
                        'amplitude criterion is skipped. Set QC.Amplitude = [] to silence this.'], ...
                        sd.InputPath, N, p.QC.Amplitude(1), max(max_amps));
                    bad_amp(:) = false;
                end
                % Spike-train criteria on the whole parent sorting when this recording is a
                % window of it: per-window checks would drop exactly the units a treatment
                % silences or drives, so the unit set would differ between segments.
                [qc_times, qc_duration] = RecordingProcessor.qcReferenceTrains(sd, N, p.QC);
                if isempty(qc_times)
                    qc_times = unit_spike_times; qc_duration = sd.Duration;
                end
                qc_rates  = cellfun(@numel, qc_times) / qc_duration;
                bad_fr    = RecordingProcessor.checkFiringRate(qc_rates(:), p.QC.FiringRate);   % column: a row broadcasts to an N x N mask
                bad_rpv   = RecordingProcessor.checkRPV(qc_times, p.QC.RPV, p.QC.RefractoryPeriod);
                [is_axon, is_noisy] = RecordingProcessor.checkWaveform(norm_wf, p.QC.Axon, p.QC.Noise, p.QC.NoiseCutout, sd.SamplingRate);
                spike_counts = cellfun(@numel, qc_times);
                min_n = p.QC.MinSpikeCount;
                bad_spike_count = spike_counts(:) < min_n;
                good = ~(bad_power | bad_amp | bad_fr | bad_rpv | is_axon | is_noisy | bad_spike_count);
            else
                good = RecordingProcessor.goodUnitMask(p.QC.GoodUnits, N);
            end
            fprintf('RecordingProcessor QC: %d/%d units passed\n', sum(good), N);

            % Optional Kilosort / Bombcell label filters. Labels are stored per sorter
            % cluster; segment templates can be a re-indexed subset of those clusters.
            use_ks = isfield(p.QC, 'KSLabel') && p.QC.KSLabel.Enable;
            use_bc = isfield(p.QC, 'Bombcell') && p.QC.Bombcell.Enable;
            ks_label = strings(N, 1);
            bc_type  = nan(N, 1);
            if use_ks || use_bc
                cluster_ids = RecordingProcessor.templateClusterIds(sd.InputPath, N);
            end
            if use_ks
                try
                    T = RecordingProcessor.readClusterTable(sd.InputPath, 'cluster_KSLabel.tsv');
                    [hit, loc] = ismember(cluster_ids, T.cluster_id);
                    ks_label(hit) = string(T.KSLabel(loc(hit)));
                    good = good(:) & ismember(ks_label, string(p.QC.KSLabel.AcceptedLabels));
                catch ME
                    warning('RecordingProcessor:ksLabel', '%s', ME.message);
                end
            end
            if use_bc
                try
                    T = RecordingProcessor.readClusterTable(sd.InputPath, 'cluster_bc_unitType.tsv');
                    [hit, loc] = ismember(cluster_ids, T.cluster_id);
                    bc_type(hit) = RecordingProcessor.bombcellTypeCode(T.bc_unitType(loc(hit)));
                    if p.QC.Bombcell.FilterUnits
                        accepted = RecordingProcessor.bombcellTypeCode(p.QC.Bombcell.AcceptedTypes);   % names ("GOOD") or codes
                        good = good(:) & ismember(bc_type, accepted);
                    end
                catch ME
                    warning('RecordingProcessor:bombcell', '%s', ME.message);
                end
            end
            if use_ks || use_bc
                fprintf('RecordingProcessor QC after label filters: %d/%d units passed\n', sum(good), N);
            end

            % Recordings with fewer than QC.N_Units good units are dropped (no units)
            if isfield(p.QC, 'N_Units') && sum(good) < p.QC.N_Units
                warning('RecordingProcessor:tooFewUnits', ['%s: %d units passed QC, fewer than ' ...
                    'QC.N_Units = %d; the recording is dropped (no units).'], ...
                    sd.InputPath, sum(good), p.QC.N_Units);
                good(:) = false;
            end

            % Build UnitData array for passing units
            good_idx = find(good);
            rec_id   = sd.RecordingID;
            sorting_id = sd.SortingID;   % processors saved before SortingID existed: resolve now
            if isempty(sorting_id) || strlength(sorting_id) == 0
                sorting_id = SpikeData.resolveSortingID(sd.InputPath, sd.ParentPath, sd.RecordingID);
            end
            ud = UnitData.empty(1, 0);   % ud(0) = ... errors when no unit passes QC
            for u = 1:numel(good_idx)
                tid = good_idx(u);
                ud(u) = UnitData(tid, ref_elec(tid), unit_spike_times{tid}, ...
                                 norm_wf(:, tid), rec_id, sd.SamplingRate, sd.Duration);
                ud(u).UnitID = RecordingProcessor.computeStableID(sd.Metadata, tid, ref_elec(tid), sorting_id);
                if use_ks, ud(u).KSLabel = ks_label(tid); end
                if use_bc, ud(u).BombcellType = bc_type(tid); end
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
                if nargin > 1
                    warning('RecordingProcessor:alreadyComputed', ['Unit features were already computed; ' ...
                        'the requested groups are ignored. Set proc.Status.UnitFeatures = "pending" to recompute.']);
                end
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
                feat_groups = ["ActivityFeatures","WaveformFeatures","RegularityFeatures","Catch22","ACG"];
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
            act_params.RefractoryPeriod = 0.005;     % LvR constant R (Shinomoto et al. 2009), not the QC window
            if isfield(proc.Parameters, 'Activity') && isfield(proc.Parameters.Activity, 'LvRRefractory')
                act_params.RefractoryPeriod = proc.Parameters.Activity.LvRRefractory;
            end
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

            % Autocorrelograms of this recording (ACG1..N). Same grid and peak
            % normalisation as Parent_ACG1..N, so FeatureStore can substitute the
            % parent version column by column (ParentFeatures = "ACG").
            if any(feat_groups == "ACG")
                acg_p = RecordingProcessor.acgParameters(proc.Parameters);
                n_bins = round(acg_p.Lag / acg_p.BinSize) * 2 + 1;
                acgs = zeros(n_bins, n_units);
                for u = 1:n_units
                    a = CellTypeClassifier.computeACG(proc.Units(u).SpikeTimes(:), ...
                        acg_p.BinSize, acg_p.Lag, sr);
                    if max(a) > 0, a = a / max(a); end
                    acgs(:, u) = a;
                end
                feature_tables{end+1} = array2table(acgs', "VariableNames", "ACG" + (1:n_bins));
            end

            new_tbl = [feature_tables{:}];
            old = proc.UnitFeatureTable;
            if width(old) > 0 && height(old) == height(new_tbl)
                % A re-run replaces this step's columns and keeps those of the other steps
                % (parent ACGs, burst participation, graph, cell-type and spatial columns),
                % whose Status stays "done"
                names = string(old.Properties.VariableNames);
                old(:, RecordingProcessor.isUnitStepColumn(names) | ...
                    ismember(names, string(new_tbl.Properties.VariableNames))) = [];
                new_tbl = [new_tbl, old];
            end
            proc.UnitFeatureTable = new_tbl;
            proc.Status.UnitFeatures = "done";
            proc.FeatureVersion.UnitFeatures = RecordingProcessor.FeatureDefinitionVersion;
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
                if nargin > 1
                    warning('RecordingProcessor:alreadyComputed', ['Parent features were already computed; ' ...
                        'the ACG parameters are ignored. Set proc.Status.ParentFeatures = "pending" to recompute.']);
                end
                return
            end
            if proc.Status.QC ~= "done"
                proc.runQC();
            end

            % ACG grid: Parameters.ACG (same as the ACG1..N columns), overridable per call
            p = parseStructParameters(RecordingProcessor.acgParameters(proc.Parameters), acg_params);

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
                % Provenance, so consumers (CellTypeClassifier) can adapt or recompute.
                proc.Units(u).FullACGBinSize = p.BinSize;
                proc.Units(u).FullACGLag     = p.Lag;
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
            proc.FeatureVersion.ParentFeatures = RecordingProcessor.FeatureDefinitionVersion;
        end

        function computeNetworkFeatures(proc)
        % COMPUTENETWORKFEATURES  Compute burst stats, regularity, Catch22 at network level.
            if proc.Status.NetworkFeatures == "done"
                return
            end

            if proc.Status.QC ~= "done"
                proc.runQC();
            end
            % Network analyses use QC-passing units only; spike_units are indices into proc.Units.
            [spike_times, spike_units] = RecordingProcessor.qcSpikes(proc);
            duration    = proc.SpikeData.Duration;
            net_struct  = struct();
            n_failed = 0;

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
                    n_failed = n_failed + 1;
                    warning('RecordingProcessor:networkFeatures', 'Regularity failed: %s', ME.message);
                end
            end

            if proc.Parameters.Analyses.Catch22
                try
                    net_struct.Catch22 = computeCatch22(spike_times, duration, proc.Parameters.Catch22.BinSize);
                catch ME
                    n_failed = n_failed + 1;
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
                    if isfield(proc.Parameters, 'Bursts') && isfield(proc.Parameters.Bursts, 'TonicMethod')
                        burst_params.TonicMethod = proc.Parameters.Bursts.TonicMethod;
                    end
                    if isfield(proc.Parameters, 'Bursts') && isfield(proc.Parameters.Bursts, 'ShapeLevel')
                        burst_params.ShapeLevel = proc.Parameters.Bursts.ShapeLevel;
                    end
                    if isfield(proc.Parameters, 'Bursts') && isfield(proc.Parameters.Bursts, 'TonicMinRate')
                        burst_params.TonicMinRate = proc.Parameters.Bursts.TonicMinRate;
                    end
                    burst_params.Outlier = proc.Parameters.Outlier;
                    burst_params.Binning = proc.Parameters.Bursts.Binning;

                    fr  = proc.UnitFeatureTable.FiringRate;
                    cv  = proc.UnitFeatureTable.CVInterSpikeInterval;
                    [bursts, burst_features] = computeBursts(spike_times, ...
                        spike_units, fr, cv, ...
                        numel(proc.Units), duration, burst_params);
                    proc.Bursts = bursts;
                    net_struct.BurstFeatures = burst_features;

                    if ~isempty(bursts.T_start)
                        n_units = numel(proc.Units);
                        edges = reshape([bursts.T_start(:)'; bursts.T_end(:)'], 1, []);
                        if issorted(edges)
                            % disjoint bursts: odd histogram bins lie inside [T_start, T_end]
                            [~, ~, bin] = histcounts(spike_times, edges);
                            in_burst = mod(bin, 2) == 1 | ismember(spike_times, bursts.T_end);
                        else
                            in_burst = false(size(spike_times));
                            for b = 1:length(bursts.T_start)
                                in_burst = in_burst | (spike_times >= bursts.T_start(b) & spike_times <= bursts.T_end(b));
                            end
                        end
                        n_spk = accumarray(spike_units(:), 1, [n_units 1]);
                        n_in  = accumarray(spike_units(:), double(in_burst(:)), [n_units 1]);
                        burst_frac = n_in ./ n_spk;
                        burst_frac(n_spk == 0) = NaN;
                        proc.UnitFeatureTable.BurstFractionSpikes = burst_frac;
                    end
                catch ME
                    n_failed = n_failed + 1;
                    warning('RecordingProcessor:networkFeatures', 'Bursts failed: %s', ME.message);
                end
            end

            % Flatten struct to a one-row table; on a re-run (after a failure) the columns of
            % later steps (graph, cell-type, spatial) are kept
            if ~isempty(fieldnames(net_struct))
                proc.NetworkFeatureTable = RecordingProcessor.mergeColumns(proc.NetworkFeatureTable, ...
                    RecordingProcessor.structToOneRowTable(net_struct));
            elseif isempty(proc.NetworkFeatureTable)
                proc.NetworkFeatureTable = table();
            end
            proc.Status.NetworkFeatures = RecordingProcessor.stepStatus(n_failed);
            proc.FeatureVersion.NetworkFeatures = RecordingProcessor.FeatureDefinitionVersion;
        end

        function computeConnectivity(proc, methods)
        % COMPUTECONNECTIVITY  Compute CCG / STTC / DDC connectivity via standalone functions.
            arguments
                proc
                methods string = proc.Parameters.Analyses.Connectivity
            end
            if isempty(methods)
                proc.Status.Connectivity = "done";
                return
            end
            if proc.Status.Connectivity == "done"
                % Methods requested now that an earlier call did not compute (e.g. "DDC" after
                % runAll) are added; the graph features are then recomputed over all of them.
                if isempty(proc.Connectivity), have = false(size(methods));
                else, have = isfield(proc.Connectivity, cellstr(methods)); end
                methods = methods(~have);
                if isempty(methods)
                    return
                end
            end
            if proc.Status.QC ~= "done"
                proc.runQC();
            end

            [st, su] = RecordingProcessor.qcSpikes(proc);   % QC units only, su = index into proc.Units (matrices are n x n)
            sr  = proc.SpikeData.SamplingRate;
            n   = numel(proc.Units);
            dur = proc.SpikeData.Duration;

            n_failed = 0;
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
                    n_failed = n_failed + 1;
                    warning('RecordingProcessor:connectivity', '%s: %s', methods(a), ME.message);
                end
            end

            % Graph features
            if ~isempty(proc.Connectivity) && ~isempty(fieldnames(proc.Connectivity))
                try
                    gp = struct();
                    if isfield(proc.Parameters, 'GraphFeatures')
                        gp = proc.Parameters.GraphFeatures;   % SmallWorldnessIterations, Types
                    end
                    [nw_graph, unit_graph] = computeGraphFeatures(proc.Connectivity, gp);
                    proc.NetworkFeatureTable = RecordingProcessor.mergeColumns(proc.NetworkFeatureTable, nw_graph);
                    if ~isempty(unit_graph) && ~isempty(proc.UnitFeatureTable)
                        proc.UnitFeatureTable = RecordingProcessor.mergeColumns(proc.UnitFeatureTable, unit_graph);
                    end
                catch ME
                    n_failed = n_failed + 1;
                    warning('RecordingProcessor:graphFeatures', '%s', ME.message);
                end
            end
            proc.Status.Connectivity = RecordingProcessor.stepStatus(n_failed);
            proc.FeatureVersion.Connectivity = RecordingProcessor.FeatureDefinitionVersion;
        end

        function computeCellTypeFeatures(proc)
        % COMPUTECELLTYPEFEATURES  Compute E/I-stratified graph and balance features.
        %   Requires CellTypeLabels to be set and Connectivity to be computed.
        %   Produces columns with _EE, _II, _EI, _IE suffixes in feature tables.
            if proc.Status.CellTypeFeatures == "done"
                return
            end
            if isempty(proc.CellTypeLabels) || all(isnan(proc.CellTypeLabels))
                proc.Status.CellTypeFeatures = "done";
                return
            end
            if isempty(proc.Connectivity) || isempty(fieldnames(proc.Connectivity))
                warning('RecordingProcessor:cellTypeFeatures', ...
                    'No connectivity data — run computeConnectivity() first.');
                proc.Status.CellTypeFeatures = "done";
                return
            end

            labels   = proc.CellTypeLabels;
            n_failed = 0;
            exc_mask = labels == 1;
            inh_mask = labels == 2;

            % Remove any stale cell-type columns before appending
            proc.NetworkFeatureTable = removeCellTypeCols(proc.NetworkFeatureTable);
            proc.UnitFeatureTable    = removeCellTypeCols(proc.UnitFeatureTable);

            % Cell-type-stratified graph features
            try
                gp = struct();
                if isfield(proc.Parameters, 'GraphFeatures') && ...
                        isfield(proc.Parameters.GraphFeatures, 'SmallWorldnessIterations')
                    gp.SmallWorldnessIterations = proc.Parameters.GraphFeatures.SmallWorldnessIterations;
                end
                [nw_ct, unit_ct] = computeCellTypeGraphFeatures(proc.Connectivity, labels, gp);
                if ~isempty(nw_ct) && ~isempty(proc.NetworkFeatureTable)
                    proc.NetworkFeatureTable = [proc.NetworkFeatureTable, nw_ct];
                elseif ~isempty(nw_ct)
                    proc.NetworkFeatureTable = nw_ct;
                end
                if ~isempty(unit_ct) && ~isempty(proc.UnitFeatureTable)
                    proc.UnitFeatureTable = [proc.UnitFeatureTable, unit_ct];
                end
            catch ME
                n_failed = n_failed + 1;
                warning('RecordingProcessor:cellTypeFeatures', ...
                    'Cell-type graph features failed: %s', ME.message);
            end

            % E/I balance features (network level)
            try
                n_exc = sum(exc_mask);
                n_inh = sum(inh_mask);
                balance = struct();
                balance.ExcitatoryFraction = n_exc / max(n_exc + n_inh, 1);
                if ~isempty(proc.UnitFeatureTable) && ...
                        ismember('FiringRate', proc.UnitFeatureTable.Properties.VariableNames)
                    fr = proc.UnitFeatureTable.FiringRate;
                    balance.MeanFiringRate_E   = mean(fr(exc_mask), 'omitnan');
                    balance.MeanFiringRate_I   = mean(fr(inh_mask), 'omitnan');
                    balance.FiringRateRatio_EI = balance.MeanFiringRate_E / ...
                        max(balance.MeanFiringRate_I, eps);
                else
                    balance.MeanFiringRate_E   = NaN;
                    balance.MeanFiringRate_I   = NaN;
                    balance.FiringRateRatio_EI = NaN;
                end
                bal_tbl = RecordingProcessor.structToOneRowTable(balance);
                if ~isempty(proc.NetworkFeatureTable)
                    proc.NetworkFeatureTable = [proc.NetworkFeatureTable, bal_tbl];
                else
                    proc.NetworkFeatureTable = bal_tbl;
                end
            catch ME
                n_failed = n_failed + 1;
                warning('RecordingProcessor:cellTypeFeatures', ...
                    'E/I balance features failed: %s', ME.message);
            end

            % Cell-type activity summaries (network-level means split by type)
            try
                feature_map = struct( ...
                    'CV2',           'CV2InterSpikeInterval', ...
                    'FanoFactor',    'FanoFactor', ...
                    'LvR',           'RevisedLocalVariation', ...
                    'ISIBimodality', 'ISIBimodality');
                summary_fields = fieldnames(feature_map);
                summary_struct = struct();
                if ~isempty(proc.UnitFeatureTable)
                    for k = 1:numel(summary_fields)
                        disp_name = summary_fields{k};
                        col_name  = feature_map.(disp_name);
                        if ismember(col_name, proc.UnitFeatureTable.Properties.VariableNames)
                            vals = proc.UnitFeatureTable.(col_name);
                        else
                            vals = NaN(numel(labels), 1);
                        end
                        summary_struct.(char("Mean" + disp_name + "_E")) = mean(vals(exc_mask), 'omitnan');
                        summary_struct.(char("Mean" + disp_name + "_I")) = mean(vals(inh_mask), 'omitnan');
                    end
                    act_tbl = RecordingProcessor.structToOneRowTable(summary_struct);
                    proc.NetworkFeatureTable = [proc.NetworkFeatureTable, act_tbl];
                end
            catch ME
                n_failed = n_failed + 1;
                warning('RecordingProcessor:cellTypeFeatures', ...
                    'Cell-type activity summaries failed: %s', ME.message);
            end

            % Cell-type burst features: E / I activity around network bursts, timed by burst
            % onset and offset, each against random "I" sets matched in firing rate
            % (computeCellTypeBurstEI). Bursts of the "ISINJoint" method carry onsets; for the
            % other methods they are detected here with detectNetworkBursts.
            try
                [st, su] = RecordingProcessor.qcSpikes(proc);
                duration = proc.SpikeData.Duration;
                bursts = proc.Bursts;
                if isempty(bursts) || ~isfield(bursts, 'Onset')
                    pb = proc.Parameters.Bursts;
                    uft = proc.UnitFeatureTable;
                    fr = NaN(numel(proc.Units), 1); cv = fr;
                    if ismember('FiringRate', uft.Properties.VariableNames), fr = uft.FiringRate; end
                    if ismember('CVInterSpikeInterval', uft.Properties.VariableNames), cv = uft.CVInterSpikeInterval; end
                    bursts = detectNetworkBursts(st, su, duration, TonicMethod=pb.TonicMethod, ...
                        TonicMinRate=pb.TonicMinRate, FiringRates=fr, CVISI=cv, ShapeLevel=pb.ShapeLevel);
                end
                cb = proc.Parameters.CellTypeBursts;
                [nw_burst, unit_burst] = computeCellTypeBurstEI(st, su, labels, bursts, duration, ...
                    NNull=cb.NNull, NullStrata=cb.NullStrata, Seed=cb.Seed, ...
                    MinUnits=cb.MinUnits, MinBursts=cb.MinBursts);
                proc.NetworkFeatureTable = [proc.NetworkFeatureTable, nw_burst];
                if ~isempty(proc.UnitFeatureTable)
                    proc.UnitFeatureTable = [proc.UnitFeatureTable, unit_burst];
                end
            catch ME
                n_failed = n_failed + 1;
                warning('RecordingProcessor:cellTypeFeatures', ...
                    'Cell-type burst features failed: %s', ME.message);
            end

            % E/I correlation structure from STTC
            try
                if isfield(proc.Connectivity, 'STTC') && isfield(proc.Connectivity.STTC, 'wu')
                    [nw_corr, unit_corr] = computeCellTypeCorrelationFeatures( ...
                        proc.Connectivity.STTC.wu, labels);
                    if ~isempty(nw_corr)
                        proc.NetworkFeatureTable = [proc.NetworkFeatureTable, nw_corr];
                    end
                    if ~isempty(unit_corr) && ~isempty(proc.UnitFeatureTable)
                        proc.UnitFeatureTable = [proc.UnitFeatureTable, unit_corr];
                    end
                end
            catch ME
                n_failed = n_failed + 1;
                warning('RecordingProcessor:cellTypeFeatures', ...
                    'E/I correlation features failed: %s', ME.message);
            end

            proc.Status.CellTypeFeatures = RecordingProcessor.stepStatus(n_failed);
            proc.FeatureVersion.CellTypeFeatures = RecordingProcessor.FeatureDefinitionVersion;
        end

        function computeSpatialAnalysis(proc)
        % COMPUTESPATIALANALYSIS  Spatial distribution features for units and cell types.
        %   Requires runQC() and computeUnitFeatures() to have run.
        %   Cell-type stratified features are computed when CellTypeLabels are available.
        %   Connectivity-based features (SpatialDecayTau) require Connectivity.STTC.
        %   Burst propagation features require Bursts to be populated.
            if proc.Status.SpatialFeatures == "done"
                return
            end
            if isempty(proc.Units)
                proc.Status.SpatialFeatures = "done";
                return
            end
            if isempty(proc.SpikeData) || isempty(proc.SpikeData.ElectrodeCoordinates)
                warning('RecordingProcessor:spatialAnalysis', ...
                    'No electrode coordinates available — skipping spatial analysis.');
                proc.Status.SpatialFeatures = "done";
                return
            end

            ec     = proc.SpikeData.ElectrodeCoordinates;
            n_failed = 0;
            labels = proc.CellTypeLabels;
            sp     = proc.Parameters;

            % Remove stale spatial columns before recomputing
            proc.NetworkFeatureTable = removeSpatialCols(proc.NetworkFeatureTable);
            proc.UnitFeatureTable    = removeSpatialCols(proc.UnitFeatureTable);

            sp_params = struct();
            if isfield(sp, 'SpatialFeatures')
                sp_params = sp.SpatialFeatures;
            end

            % --- Core spatial features (A + B + E + optional G) ---
            try
                [nw_sp, unit_sp] = computeSpatialFeatures( ...
                    ec, proc.Units, labels, proc.UnitFeatureTable, sp_params);
                RecordingProcessor.appendFeatures(proc, nw_sp, unit_sp);
            catch ME
                n_failed = n_failed + 1;
                warning('RecordingProcessor:spatialAnalysis', ...
                    'Core spatial features failed: %s', ME.message);
            end

            % --- Distance-dependent FC decay (D) ---
            try
                if ~isempty(proc.Connectivity) && ~isempty(fieldnames(proc.Connectivity))
                    [nw_conn, unit_conn] = computeSpatialConnectivityFeatures( ...
                        ec, proc.Units, proc.Connectivity, labels, sp_params);
                    RecordingProcessor.appendFeatures(proc, nw_conn, unit_conn);
                end
            catch ME
                n_failed = n_failed + 1;
                warning('RecordingProcessor:spatialAnalysis', ...
                    'Spatial connectivity features failed: %s', ME.message);
            end

            % --- Spatial E/I balance map (F) ---
            try
                if ~isempty(labels) && any(~isnan(labels))
                    [nw_ei, ~] = computeSpatialEIBalance( ...
                        ec, proc.Units, labels, sp_params);
                    RecordingProcessor.appendFeatures(proc, nw_ei, table());
                end
            catch ME
                n_failed = n_failed + 1;
                warning('RecordingProcessor:spatialAnalysis', ...
                    'Spatial E/I balance failed: %s', ME.message);
            end

            % --- Burst propagation (H) ---
            try
                if ~isempty(proc.Bursts) && isfield(proc.Bursts, 'T_start') ...
                        && ~isempty(proc.Bursts.T_start)
                    [nw_burst, ~] = computeSpatialBurstFeatures( ...
                        ec, proc.Units, proc.Bursts, ...
                        proc.SpikeData.SpikeTimes, proc.SpikeData.SpikeUnits, ...
                        labels, sp_params);
                    RecordingProcessor.appendFeatures(proc, nw_burst, table());
                end
            catch ME
                n_failed = n_failed + 1;
                warning('RecordingProcessor:spatialAnalysis', ...
                    'Spatial burst features failed: %s', ME.message);
            end

            proc.Status.SpatialFeatures = RecordingProcessor.stepStatus(n_failed);
            proc.FeatureVersion.SpatialFeatures = RecordingProcessor.FeatureDefinitionVersion;
        end

        function runAll(proc)
        % RUNALL  Convenience: run all enabled analyses in sequence.
            proc.runQC();
            if ~isempty(proc.Units)
                proc.computeUnitFeatures();
                proc.computeParentFeatures();
                proc.computeNetworkFeatures();
                proc.computeConnectivity();
                proc.computeCellTypeFeatures();
                proc.computeSpatialAnalysis();
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
            CellTypeLabels        = proc.CellTypeLabels;      %#ok<PROP>
            Parameters            = proc.Parameters;          %#ok<PROP>
            Status                = proc.Status;              %#ok<PROP>
            FeatureVersion        = proc.FeatureVersion;      %#ok<PROP>
            builtin('save', file_path, 'SpikeDataStruct', 'UnitsStructArray', ...
                'UnitFeatureTable', 'NetworkFeatureTable', ...
                'Connectivity', 'Bursts', 'CellTypeLabels', 'Parameters', 'Status', 'FeatureVersion');
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
            if isfield(s, 'CellTypeLabels')
                proc.CellTypeLabels = s.CellTypeLabels;
            end
            if isfield(s, 'FeatureVersion')
                proc.FeatureVersion = s.FeatureVersion;
            end
            if isfield(s, 'Parameters')
                % Parameters added since the file was saved take their defaults (steps read
                % them unguarded); stored values, including retired fields, are kept as saved.
                proc.Parameters = RecordingProcessor.fillMissingParameters( ...
                    RecordingProcessor.returnDefaultParams(), s.Parameters);
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
                if ~isfield(proc.Status, 'CellTypeFeatures')
                    proc.Status.CellTypeFeatures = "pending";
                end
                if ~isfield(proc.Status, 'SpatialFeatures')
                    proc.Status.SpatialFeatures = "pending";
                end
                % Upgrade legacy boolean Status fields to tri-state strings.
                proc.Status = RecordingProcessor.upgradeLegacyStatus(proc.Status);
            else
                proc.Status = RecordingProcessor.emptyStatus();
            end
        end

        function procs = loadMany(file_paths)
        % LOADMANY  Load multiple processors; in parallel only if a pool is already open.
        %   Starting a pool just for loading costs more than it saves, and concurrent MATLAB
        %   jobs that each auto-start a pool collide in the shared local cluster folder.
            if ischar(file_paths) || isstring(file_paths)
                file_paths = cellstr(file_paths);
            end
            N = numel(file_paths);
            procs(N) = RecordingProcessor();
            n_workers = 0;
            if exist('gcp', 'file') == 2
                pool = gcp('nocreate');
                if ~isempty(pool), n_workers = pool.NumWorkers; end
            end
            parfor (i = 1:N, n_workers)
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
                % After QC, MEArecording.updateSpikeTimes stores the spike owners as compact
                % unit indices (1..numel(Units)); the new pipeline identifies spikes by
                % TemplateID (qcSpikes, burst and spatial features), so translate them.
                su = proc.SpikeData.SpikeUnits;
                if ~isempty(su) && all(su >= 1 & su <= numel(proc.Units) & su == round(su))
                    tids = [proc.Units.TemplateID];
                    proc.SpikeData.SpikeUnits = reshape(tids(su), size(su));
                end
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
            params.Analyses.SingleCell   = true;
            params.Analyses.Regularity   = true;
            params.Analyses.Bursts       = true;
            params.Analyses.Connectivity = ["CCG","STTC"];
            % catch22 needs the compiled feature MEX files (run mexAll in
            % Toolboxes/catch22/wrap_Matlab); default on only if they are present.
            % catch22_all itself is an .m wrapper, so test one compiled feature.
            params.Analyses.Catch22      = exist('catch22_DN_HistogramMode_5', 'file') == 3;

            params.Outlier.Method          = "median";
            params.Outlier.ThresholdFactor = 5;

            params.QC.LSB              = 6.2;
            params.QC.Amplitude        = [20 1000];
            params.QC.FiringRate       = [0.01 10];
            params.QC.RPV              = 0.02;
            params.QC.RefractoryPeriod = 0.002;
            params.QC.Axon             = 0.8;
            params.QC.Noise            = 14;
            params.QC.NoiseCutout      = [-1 1];
            params.QC.PowerCutoff      = 1.3;
            params.QC.N_Units          = 10;
            params.QC.MinSpikeCount    = 10;
            % "parent": firing rate, spike count and RPV of a window of a longer sorting
            % (split_sortings / segment_N / ParentPath) are checked on the whole sorting, so
            % every segment keeps the same units; "recording": on this window only
            params.QC.RateReference    = "parent";
            params.QC.GoodUnits        = [];
            params.QC.KSLabel.Enable         = false;
            params.QC.KSLabel.AcceptedLabels = "good";
            params.QC.Bombcell.Enable        = false;
            params.QC.Bombcell.Path          = '';
            params.QC.Bombcell.FilterUnits   = true;
            params.QC.Bombcell.AcceptedTypes = [1, 2];
            params.QC.Bombcell.Overrides     = struct();

            params.Activity.FanoBinWidth = 0.1;
            params.Activity.LvRRefractory = 0.005;   % s, refractory constant R of the revised local variation (Shinomoto et al. 2009)

            params.Bursts.Method      = "ISIN";      % legacy; "ISINJoint" = detectNetworkBursts, "logISI"
            params.Bursts.TonicMethod = "gini";      % "ISINJoint": unit removal before detection (with TonicMinRate)
            params.Bursts.ShapeLevel  = 0.25;        % burst onset / offset at this fraction of the peak (burstShape)
            params.Bursts.TonicMinRate = "median";   % "ISINJoint": units below this rate (Hz or "median") are never removed
            params.Bursts.MergeFactor = 0.5;
            params.Bursts.Binning     = 0.1;

            % Cell-type burst features (computeCellTypeBurstEI): rate-matched null draws and
            % the minimum bursts / units per type below which they are NaN
            params.CellTypeBursts.NNull      = 200;
            params.CellTypeBursts.NullStrata = 5;
            params.CellTypeBursts.Seed       = 1;
            params.CellTypeBursts.MinUnits   = 3;
            params.CellTypeBursts.MinBursts  = 5;

            params.Regularity.Binning = 0.1;
            params.Regularity.N_peaks = 20;

            % Unit autocorrelograms (ACG1..N, and Parent_ACG1..N from the parent sorting):
            % same grid as CellTypeClassifier.Harmonization, peak-normalised.
            params.ACG.BinSize = 0.0005;
            params.ACG.Lag     = 0.1;

            params.CCG.BinSize  = 0.001;
            params.CCG.Duration = 0.1;
            params.CCG.Conv_w   = 0.010 / params.CCG.BinSize;
            params.CCG.Alpha    = 0.001;
            params.CCG.InhibitoryCriterion = "strict";   % or "cellexplorer" (CellExplorer / previous behaviour; see computeConnectivityCCG)

            params.DDC.BinSize     = 0.001;
            params.DDC.Threshold   = 0;
            params.DDC.N_Surrogates = 100;
            params.DDC.Percentile  = 95;

            params.STTC.MaxLag          = 0.05;
            params.STTC.SurrogateMethod = 'jitter';
            params.STTC.JitterSD        = [];      % s, surrogate jitter SD; empty = 2 x MaxLag (must exceed MaxLag)
            params.STTC.N_Surrogates    = 100;
            params.STTC.EdgeTest        = "fdr";   % "fdr": z vs surrogates + BH across pairs; "percentile": per-pair
            params.STTC.FDR             = 0.05;
            params.STTC.Percentile      = 95;      % EdgeTest "percentile" only

            params.GraphFeatures.SmallWorldnessIterations = 100;
            % Types: "binary" (default), "weighted", or ["binary","weighted"]. Weighted
            % graphs use the link weights of the significant links (DDC: |coefficient|,
            % STTC: tiling coefficient); CCG has binary links only.
            params.GraphFeatures.Types = "binary";

            params.Catch22.BinSize = 0.1;

            params.SpatialFeatures.NNk                       = 5;
            params.SpatialFeatures.RipleyEnable              = false;
            params.SpatialFeatures.RipleyMaxRadius           = 500;
            params.SpatialFeatures.RipleyNBins               = 50;
            params.SpatialFeatures.RipleyNSurrogates         = 100;
            params.SpatialFeatures.KernelBandwidth           = 100;
            params.SpatialFeatures.BurstMinUnitsForWavefront = 5;

            params.Save.Flag      = false;
            params.Save.Path      = [];
            params.Save.Overwrite = false;
        end

        function applyLabelsFromClassifier(proc_array, ctc)
        % APPLYLABELSFROMCLASSIFIER  Copy UnitLabels from a CellTypeClassifier to processors.
        %   Maps labels by UnitID within the same culture (ctc.Parameters.CultureKeys): a unit
        %   keeps its label in every recording of its culture (e.g. all doses), while UnitIDs
        %   that collide across cultures/sortings cannot pick up another unit's label.
        %   Resets CellTypeFeatures and SpatialFeatures to "pending" (both depend on labels).
            arguments
                proc_array RecordingProcessor
                ctc        CellTypeClassifier
            end
            ut           = ctc.FeatureStore.UnitTable;
            keys_fn      = ctc.Parameters.CultureKeys;
            ctc_cultures = FeatureStore.getCultureIDsForUnits(ut.RecordingID, ctc.FeatureStore.MetadataTable, keys_fn);
            ctc_keys     = string(ut.UnitID) + "|" + ctc_cultures;
            ctc_labels   = ctc.UnitLabels;
            for i = 1:numel(proc_array)
                p = proc_array(i);
                if isempty(p.Units), continue; end
                proc_culture = FeatureStore.buildCultureIDs(struct2table(p.SpikeData.Metadata, 'AsArray', true), keys_fn);
                proc_keys    = string({p.Units.UnitID}) + "|" + proc_culture;
                labels       = NaN(1, numel(p.Units));
                [found, loc] = ismember(proc_keys, ctc_keys);
                labels(found) = ctc_labels(loc(found));
                p.CellTypeLabels = labels;
                p.Status.CellTypeFeatures = "pending";
                p.Status.SpatialFeatures  = "pending";
            end
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Static, Access = private)

        function p = fillMissingParameters(defaults, p)
        % Add fields of defaults that p lacks, at any depth; existing fields of p are kept.
            for f = string(fieldnames(defaults))'
                if ~isfield(p, f)
                    p.(f) = defaults.(f);
                elseif isstruct(defaults.(f)) && isscalar(defaults.(f)) && ...
                        isstruct(p.(f)) && isscalar(p.(f))
                    p.(f) = RecordingProcessor.fillMissingParameters(defaults.(f), p.(f));
                end
            end
        end

        function [st, su] = qcSpikes(proc)
        % QCSPIKES  Spike times of QC-passing units, with unit labels 1..numel(proc.Units).
        %   SpikeData holds every template; network and connectivity analyses must not see
        %   rejected templates, and their matrices are indexed by position in proc.Units.
            [keep, su] = ismember(proc.SpikeData.SpikeUnits, [proc.Units.TemplateID]);
            st = proc.SpikeData.SpikeTimes(keep);
            su = su(keep);
        end

        function c = splitByUnit(st, su, N)
        % SPLITBYUNIT  {1 x N} spike times per unit label 1..N (one sort instead of N scans).
            st = st(:); su = su(:);
            keep = su >= 1 & su <= N & su == round(su);
            if N == 0, c = cell(1, 0); return, end
            [s, o] = sort(su(keep));            % stable: spike order within a unit is kept
            t = st(keep);
            c = mat2cell(t(o), accumarray(s, 1, [N 1]), 1)';
        end

        function good = goodUnitMask(gu, N)
        % GOODUNITMASK  QC.GoodUnits as an (N x 1) logical mask: a logical (or 0/1) mask over the
        %   N templates, or 1-based template indices.
            gu = gu(:);
            if islogical(gu) || (numel(gu) == N && all(gu == 0 | gu == 1) && any(gu == 0))
                assert(numel(gu) == N, 'RecordingProcessor:goodUnits', ...
                    'QC.GoodUnits mask has %d entries for %d templates.', numel(gu), N);
                good = logical(gu);
            else
                assert(all(gu >= 1 & gu <= N & gu == round(gu)), 'RecordingProcessor:goodUnits', ...
                    'QC.GoodUnits indices must be template numbers 1..%d.', N);
                good = false(N, 1);
                good(gu) = true;
            end
        end

        function tf = isUnitStepColumn(names)
        % ISUNITSTEPCOLUMN  Columns written by computeUnitFeatures (any feature group).
            names = string(names);
            wf = ["AUC_peak_1","AUC_trough","AUC_peak_2","Rise","Decay","HalfWidth", ...
                  "Asymmetry","T2Pdelay","T2Pratio"];
            tf = ismember(names, [unitFeatureNames("act"), unitFeatureNames("reg"), wf]) | ...
                 startsWith(names, "SC_") | ~cellfun(@isempty, regexp(cellstr(names), '^(Waveform|ACG)\d+$', 'once'));
        end

        function T = mergeColumns(T, new)
        % MERGECOLUMNS  Append the columns of new to T, replacing columns of the same name
        %   (a re-run step overwrites its own columns and keeps those of the other steps).
            if isempty(new) || width(new) == 0
                return
            end
            if isempty(T) || width(T) == 0
                T = new;
                return
            end
            T(:, ismember(T.Properties.VariableNames, new.Properties.VariableNames)) = [];
            T = [T, new];
        end

        function s = stepStatus(n_failed)
        % STEPSTATUS  "done", or "failed" when a sub-analysis of the step raised an error (the
        %   step is then re-run by the next call and FeatureStore reports it).
            if n_failed > 0
                s = "failed";
            else
                s = "done";
            end
        end

        function s = emptyStatus()
            s.QC               = "pending";
            s.UnitFeatures     = "pending";
            s.ParentFeatures   = "pending";
            s.NetworkFeatures  = "pending";
            s.Connectivity     = "pending";
            s.CellTypeFeatures = "pending";
            s.SpatialFeatures  = "pending";
        end

        function proc = appendFeatures(proc, nw_tbl, unit_tbl)
        % Append non-empty feature tables to existing tables.
            if ~isempty(nw_tbl)
                if ~isempty(proc.NetworkFeatureTable)
                    proc.NetworkFeatureTable = [proc.NetworkFeatureTable, nw_tbl];
                else
                    proc.NetworkFeatureTable = nw_tbl;
                end
            end
            if ~isempty(unit_tbl) && ~isempty(proc.UnitFeatureTable)
                proc.UnitFeatureTable = [proc.UnitFeatureTable, unit_tbl];
            end
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

        function uid = computeStableID(metadata, template_id, ref_electrode, sorting_id)
        % Compute deterministic StableID from recording metadata + unit identity.
        %   The sorting identity (SpikeData.SortingID) keeps units of separately sorted
        %   recordings of the same culture apart (template k of one day's sorting is unrelated
        %   to template k of another day's), while segments of one sorting share their IDs.
            if nargin < 4, sorting_id = ''; end
            chip_id     = '';
            plating_date = '';
            if isfield(metadata, 'ChipID'),      chip_id     = char(string(metadata.ChipID));     end
            if isfield(metadata, 'PlatingDate'), plating_date = char(string(metadata.PlatingDate)); end
            uid = paramHash(struct('ChipID', chip_id, 'PlatingDate', plating_date, ...
                'SortingID', char(string(sorting_id)), ...
                'TemplateID', template_id, 'ReferenceElectrode', ref_electrode));
        end

        function T = readClusterTable(input_path, file_name)
        % READCLUSTERTABLE  Per-cluster TSV (phy/Bombcell) from the segment or its sorter folder.
            candidates = { fullfile(input_path, file_name), ...
                           fullfile(fileparts(input_path), file_name) };
            found = candidates(cellfun(@(f) exist(f, 'file') == 2, candidates));
            if isempty(found)
                error('DeePhys:noClusterTable', '%s not found in %s or its parent folder', file_name, input_path);
            end
            T = readtable(found{1}, 'FileType', 'text', 'Delimiter', '\t');
        end

        function [times, duration] = qcReferenceTrains(sd, n_templates, qc)
        % QCREFERENCETRAINS  Parent-sorting spike trains (cell, one per template) and duration
        %   for the spike-train QC criteria of a recording that is a window of a longer
        %   sorting; {} when QC.RateReference is "recording" or no parent sorting is found.
        %   Parent: the sorting named in split_sortings' split_window.txt, the sorter output
        %   holding a segment_N / split window folder, else SpikeData.ParentPath.
            times = {}; duration = NaN;
            if ~isfield(qc, 'RateReference') || string(qc.RateReference) ~= "parent"
                return
            end
            in_path = char(sd.InputPath);
            parent = "";
            marker = fullfile(in_path, 'split_window.txt');
            if isfile(marker)
                tok = regexp(fileread(marker), '(?m)^parent\s*=\s*(.+?)\s*$', 'tokens', 'once');
                if ~isempty(tok) && isfile(fullfile(tok{1}, 'spike_times.npy'))
                    parent = string(tok{1});
                end
            end
            if parent == ""
                [up, leaf] = fileparts(in_path);
                if (isfile(marker) || ~isempty(regexp(leaf, '^segment_\d+$', 'once'))) && ...
                        isfile(fullfile(up, 'spike_times.npy'))
                    parent = string(up);
                elseif strlength(sd.ParentPath) > 0 && isfile(fullfile(sd.ParentPath, 'spike_times.npy'))
                    parent = sd.ParentPath;
                end
            end
            if parent == "", return, end
            try
                [trains, duration] = ParentSpikeLoader.loadSpikeTrains(parent, sd.SamplingRate);
            catch ME
                warning('RecordingProcessor:parentQC', ['%s: parent sorting %s unreadable (%s); ' ...
                    'spike-train QC uses this window only.'], sd.InputPath, parent, ME.message);
                return
            end
            times = repmat({zeros(0, 1)}, n_templates, 1);
            n = min(n_templates, numel(trains));
            times(1:n) = trains(1:n);
            % A window holds a subset of its parent's spikes per template; more spikes in the
            % window than in the parent means the template numbering differs
            su = sd.SpikeUnits(:);
            su = su(su >= 1 & su <= n_templates & su == round(su));
            own = accumarray(su, 1, [n_templates 1]);
            if any(own > cellfun(@numel, times))
                warning('RecordingProcessor:parentQC', ['%s: template numbering does not match the ' ...
                    'parent sorting %s; spike-train QC uses this window only.'], sd.InputPath, parent);
                times = {}; duration = NaN;
                return
            end
            fprintf('RecordingProcessor QC: spike-train criteria from parent sorting %s (%.0f s)\n', ...
                parent, duration);
        end

        function ids = templateClusterIds(input_path, n_templates)
        % TEMPLATECLUSTERIDS  Sorter cluster ID of each template (1..n_templates).
        %   Segmented recordings store a re-indexed subset of the sorter's clusters; the
        %   mapping is qc_output/phy_ids.npy. Without it, template k is cluster k-1.
            ids = (0:n_templates - 1)';
            candidates = { fullfile(input_path, 'phy_ids.npy'), ...
                           fullfile(fileparts(input_path), 'qc_output', 'phy_ids.npy'), ...
                           fullfile(input_path, 'qc_output', 'phy_ids.npy') };
            found = candidates(cellfun(@isfile, candidates));
            if isempty(found), return; end
            phy = double(readNPY(found{1}));
            if numel(phy) == n_templates
                ids = phy(:);
            else
                warning('RecordingProcessor:phyIds', ...
                    '%s has %d entries but the recording has %d templates; assuming template k = cluster k-1.', ...
                    found{1}, numel(phy), n_templates);
            end
        end

        function codes = bombcellTypeCode(types)
        % BOMBCELLTYPECODE  Bombcell unit types as codes (0 noise, 1 good, 2 MUA,
        %   3 non-somatic good, 4 non-somatic MUA); numeric input is passed through.
            if isnumeric(types), codes = double(types(:)); return; end
            names = upper(strtrim(string(types(:))));
            codes = nan(size(names));
            codes(names == "NOISE")         = 0;
            codes(names == "GOOD")          = 1;
            codes(names == "MUA")           = 2;
            codes(names == "NON-SOMA" | names == "NON-SOMA GOOD") = 3;
            codes(names == "NON-SOMA MUA")  = 4;
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
                    % The window is fixed in ms, so the number of samples (and of possible
                    % slope changes) grows with the sampling rate; the threshold is defined
                    % for 10 kHz (MaxWell) and scaled accordingly.
                    is_noisy  = n_changes > noise_thresh * SAMPLES_PER_MS / 10 | bad_peak';
                end
            end
            fprintf('QC noisy: %d\n', sum(is_noisy));
        end

        function p = acgParameters(params)
        % ACG grid (BinSize, Lag in s) from Parameters.ACG; defaults for processors saved
        % before the field existed.
            p = struct('BinSize', 0.0005, 'Lag', 0.1);
            if isfield(params, 'ACG')
                p = parseStructParameters(p, params.ACG);
            end
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

function tbl = removeSpatialCols(tbl)
% Remove previously computed spatial columns so they can be recomputed.
    if isempty(tbl), return; end
    cols = string(tbl.Properties.VariableNames);
    sp_prefixes = [ ...
        "ConvexHullArea", "ChipCoverage", ...
        "MeanPairwiseDistance", "CentroidSpread", ...
        "SpatialMixingIndex", "MeanFractionExcNN_E", "MeanFractionExcNN_I", ...
        "FractionExcNN", "FractionInhNN", "DistFromCentroid", ...
        "SpatialFRMoransI", "CenterPeripheryFR_ratio", ...
        "RipleysL_max", "ClusterScale", ...
        "SpatialDecayTau", "DistanceFCCorrelation", "MeanFC_NearNeighbors", ...
        "SpatialEIVariability", "SpatialEIMoransI", ...
        "BurstOriginDispersion", "MeanBurstPropagationSpeed", "BurstOriginExcFraction"];
    remove = false(1, numel(cols));
    for nm = sp_prefixes
        remove = remove | startsWith(cols, nm);
    end
    tbl(:, remove) = [];
end

function tbl = removeCellTypeCols(tbl)
% Remove previously computed cell-type columns so they can be recomputed.
    if isempty(tbl), return; end
    cols = string(tbl.Properties.VariableNames);
    ct_patterns = ["_EE", "_II", "_EI", "_IE"];
    ct_prefixes = [ ...
        "CrossTypeDensity", "WithinTypeConnections", "CrossTypeConnections", ...
        "ExcitatoryFraction", "MeanFiringRate_E", "MeanFiringRate_I", "FiringRateRatio_EI", ...
        "MeanCV2_E", "MeanCV2_I", "MeanFanoFactor_E", "MeanFanoFactor_I", ...
        "MeanLvR_E", "MeanLvR_I", "MeanISIBimodality_E", "MeanISIBimodality_I", ...
        "BurstLeadTime_E", "BurstLeadTime_I", "BurstLeadFraction", ...
        "MeanBurstFraction_E", "MeanBurstFraction_I", ...
        "IntraBurstRate_E", "IntraBurstRate_I", ...
        "BurstParticipation_E", "BurstParticipation_I", "BurstParticipationRate", ...
        "BurstEI_", "BurstUnit", ...
        "MeanSTTC_EE", "MeanSTTC_II", "MeanSTTC_EI", "STTCSynchronyIndex", ...
        "MeanSTTC_WithinType", "MeanSTTC_CrossType"];
    remove = false(1, numel(cols));
    for p = ct_patterns
        remove = remove | endsWith(cols, p);
    end
    for p = ct_prefixes
        remove = remove | startsWith(cols, p);
    end
    tbl(:, remove) = [];
end
