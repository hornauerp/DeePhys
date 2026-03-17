function [waveforms, acgs, sr_wf] = extractUnitWaveformsAndACGs(ctc, unit_list)
% EXTRACTUNITWAVEFORMSANDACGS  Extract raw waveforms and ACGs from DeePhys Unit objects.
%
% Waveforms are taken from unit.ReferenceWaveform at the native recording SR.
%
% ACG strategy (controlled by ctc.Parameters.Harmonization.ACGSource):
%   - "FullACG" (default): uses unit.FullACG (high-res, distinct from unit.ACG)
%   - "ACG":               uses unit.ACG (standard-res, coarser params)
%   In both cases, the cached value is used directly when its bin count matches
%   the Harmonization params; otherwise units are grouped by parent MEArecording
%   and one batch CCG() call is made per recording using already-in-memory spike
%   times (no .npy re-reads).
%
% INPUTS:
%   ctc       - CellTypeClassifier (provides Harmonization parameters)
%   unit_list - (1 x N) Unit array
%
% OUTPUTS:
%   waveforms - (N_samples x N_units) reference waveforms at native recording SR
%   acgs      - (N_bins x N_units)   ACGs at Harmonization.ACGBinSize / ACGLag
%   sr_wf     - native recording sampling rate in Hz

arguments
    ctc       CellTypeClassifier
    unit_list (1,:) Unit
end

ph   = ctc.Parameters.Harmonization;
N    = numel(unit_list);

% Collect per-unit sampling rates and assert uniformity.
% If rates differ, resample each waveform to the first unit's SR before stacking.
unit_srs = arrayfun(@(u) getSamplingRate(u), unit_list);
sr_wf    = unit_srs(1);

if all(unit_srs == sr_wf)
    waveforms = [unit_list.ReferenceWaveform];   % (N_samples x N_units), fast path
else
    warning('extractUnitWaveformsAndACGs:mixedSamplingRates', ...
        'Units have mixed waveform sampling rates (%s Hz). Resampling all to %g Hz.', ...
        strjoin(string(unique(unit_srs)), ', '), sr_wf);
    all_wf = cell(1, N);
    for u = 1:N
        wf = unit_list(u).ReferenceWaveform(:);
        if unit_srs(u) ~= sr_wf
            [p, q]  = rat(sr_wf / unit_srs(u), 1e-6);
            wf      = resample(wf, p, q);
        end
        all_wf{u} = wf;
    end
    max_len   = max(cellfun(@numel, all_wf));
    waveforms = zeros(max_len, N);
    for u = 1:N
        waveforms(1:numel(all_wf{u}), u) = all_wf{u};
    end
end

% ACG extraction: use pre-computed ACG/FullACG when harmonization params match
% its dimensions; recompute per recording (batch CCG) when they differ.
acg_source  = ph.ACGSource;   % "FullACG" or "ACG"
n_bins_harm = round(2 * ph.ACGLag / ph.ACGBinSize) + 1;
acgs        = zeros(n_bins_harm, N);

try
    if acg_source == "FullACG"
        n_bins_cached = length(unit_list(1).FullACG);
    else
        n_bins_cached = length(unit_list(1).ACG);
    end
    % Check both bin count AND actual parameters (bin size + lag).
    % Two different (BinSize, Lag) pairs can yield the same bin count.
    params_match = false;
    if n_bins_cached == n_bins_harm
        rec1 = unit_list(1).MEArecording;
        if ~isempty(rec1)
            if acg_source == "FullACG" && isfield(rec1.Connectivity, 'FullCCG') ...
                    && isfield(rec1.Connectivity.FullCCG, 'args')
                cached_args = rec1.Connectivity.FullCCG.args;
                params_match = abs(cached_args.binsize - ph.ACGBinSize) < 1e-9 ...
                    && abs(cached_args.duration - 2 * ph.ACGLag) < 1e-9;
            elseif acg_source == "ACG" && isfield(rec1.Connectivity, 'CCG') ...
                    && isfield(rec1.Connectivity.CCG, 'args')
                cached_args = rec1.Connectivity.CCG.args;
                params_match = abs(cached_args.binsize - ph.ACGBinSize) < 1e-9 ...
                    && abs(cached_args.duration - 2 * ph.ACGLag) < 1e-9;
            else
                % No stored args — fall back to bin count match only
                params_match = true;
            end
        else
            % No MEArecording reference — trust bin count
            params_match = true;
        end
    end
    use_cached = params_match;
catch
    use_cached    = false;
    n_bins_cached = 0;
end

if use_cached
    for u = 1:N
        if acg_source == "FullACG"
            acgs(:, u) = unit_list(u).FullACG;
        else
            acgs(:, u) = unit_list(u).ACG;
        end
    end
else
    fprintf('Harmonization ACG params (%i bins) differ from cached %s (%i bins) — recomputing from parent recording\n', ...
        n_bins_harm, acg_source, n_bins_cached);

    % Ensure CCGHeart MEX is compiled (required by CCG).
    % Check for the binary directly — which('CCGHeart') would find a .m stub
    % and return non-empty even when the MEX is missing.
    func_dir = fileparts(mfilename('fullpath'));
    if ~isfile(fullfile(func_dir, ['CCGHeart.' mexext()]))
        prev_dir = cd(func_dir);
        cleanup  = onCleanup(@() cd(prev_dir));
        mex('CCGHeart.c');
    end

    % ── FullACG recompute: read from ParentRecordingPath (.npy) ────────────
    % The FullACG uses all spikes from the concatenated parent recording
    % (across all DIVs/segments), not just the current segment's spike times.
    % Group units by ParentRecordingPath so one CCG call covers all segments.
    if acg_source == "FullACG"
        parent_paths  = {};   % unique parent paths
        parent_groups = {};   % {i} = global unit indices sharing parent_paths{i}
        parent_tids   = {};   % {i} = TemplateIDs for those units
        parent_srs    = [];   % sampling rate per parent group

        for u = 1:N
            rec = unit_list(u).MEArecording;
            if isempty(rec)
                warning('extractUnitWaveformsAndACGs:noParent', ...
                    'Unit %d has no MEArecording reference — skipping FullACG recompute.', u);
                continue
            end
            pp = rec.ParentRecordingPath;
            found = false;
            for p = 1:numel(parent_paths)
                if strcmp(parent_paths{p}, pp)
                    parent_groups{p}(end+1) = u;
                    parent_tids{p}(end+1)   = unit_list(u).TemplateID;
                    found = true;
                    break;
                end
            end
            if ~found
                parent_paths{end+1}  = pp;
                parent_groups{end+1} = u;
                parent_tids{end+1}   = unit_list(u).TemplateID;
                parent_srs(end+1)    = rec.RecordingInfo.SamplingRate;
            end
        end

        for p = 1:numel(parent_paths)
            pp    = parent_paths{p};
            g_idx = parent_groups{p};   % global indices into unit_list
            tids  = parent_tids{p};     % TemplateIDs for units in this group
            sr_p  = parent_srs(p);

            % Read full concatenated spike data from parent Kilosort output
            npy_times = fullfile(pp, 'spike_times.npy');
            npy_temps = fullfile(pp, 'spike_templates.npy');
            if ~isfile(npy_times) || ~isfile(npy_temps)
                warning('extractUnitWaveformsAndACGs:missingNPY', ...
                    'Parent .npy files not found at %s — falling back to segment spike times.', pp);
                acgs = recomputeACGsFromSegments(unit_list, g_idx, acgs, ph, n_bins_harm);
                continue
            end

            spike_times_raw = readNPY(npy_times);
            spike_templates = readNPY(npy_temps);

            % Filter to good template IDs (unique across units in this parent)
            unique_tids = unique(tids);
            keep = ismember(spike_templates, unique_tids);
            spike_times_s = double(spike_times_raw(keep)) / sr_p;
            spike_templates_kept = spike_templates(keep);

            % Assign compact local unit IDs via findgroups on kept templates
            [unit_ids_local, group_tids] = findgroups(spike_templates_kept);

            if isempty(spike_times_s)
                continue
            end

            % Sort by time (required by CCG)
            [spike_times_s, sort_idx] = sort(spike_times_s);
            unit_ids_local = unit_ids_local(sort_idx);

            [ccg_3d, ~] = CCG(spike_times_s, unit_ids_local, ...
                'binSize', ph.ACGBinSize, ...
                'duration', 2 * ph.ACGLag, ...
                'Fs', 1 / sr_p);

            % Map CCG diagonals back to unit_list indices via TemplateID
            for k = 1:numel(g_idx)
                tid = tids(k);
                lid = find(group_tids == tid, 1);
                if ~isempty(lid)
                    acgs(:, g_idx(k)) = double(ccg_3d(:, lid, lid));
                end
            end

            fprintf('  FullACG recomputed from %s (%i units, %i spikes)\n', ...
                pp, numel(g_idx), numel(spike_times_s));
        end

    else
        % ── Regular ACG recompute: use per-segment spike times ─────────────
        % For ACGSource=="ACG", segment-level spikes are appropriate.
        acgs = recomputeACGsFromSegments(unit_list, 1:N, acgs, ph, n_bins_harm);
    end
end
end

function acgs = recomputeACGsFromSegments(unit_list, indices, acgs, ph, n_bins_harm)
% Recompute ACGs from per-segment in-memory spike times.
% Used for ACGSource=="ACG" or as fallback when parent .npy files are missing.

    % Group units by parent MEArecording for one batch CCG call per recording.
    rec_list   = {};
    rec_groups = {};
    for ii = 1:numel(indices)
        u   = indices(ii);
        rec = unit_list(u).MEArecording;
        if isempty(rec)
            continue
        end
        found = false;
        for r = 1:numel(rec_list)
            if rec_list{r} == rec
                rec_groups{r}(end+1) = u;
                found = true;
                break;
            end
        end
        if ~found
            rec_list{end+1}   = rec;
            rec_groups{end+1} = u;
        end
    end

    for r = 1:numel(rec_list)
        g_idx = rec_groups{r};
        sr_r  = unit_list(g_idx(1)).getSamplingRate();

        % Build concatenated spike train with compact local unit IDs
        all_times    = [];
        all_unit_ids = [];
        local_ids    = zeros(1, numel(g_idx));
        local_ctr    = 0;
        for k = 1:numel(g_idx)
            st = unit_list(g_idx(k)).SpikeTimes;
            if ~isempty(st)
                local_ctr      = local_ctr + 1;
                local_ids(k)   = local_ctr;
                all_times      = [all_times;    st(:)];
                all_unit_ids   = [all_unit_ids; local_ctr * ones(numel(st), 1)];
            end
        end

        if local_ctr == 0
            continue
        end

        [~, sort_idx]  = sort(all_times);
        all_times      = all_times(sort_idx);
        all_unit_ids   = all_unit_ids(sort_idx);

        [ccg_3d, ~] = CCG(all_times, all_unit_ids, ...
            'binSize', ph.ACGBinSize, ...
            'duration', 2 * ph.ACGLag, ...
            'Fs', 1 / sr_r);

        for k = 1:numel(g_idx)
            if local_ids(k) > 0
                lid = local_ids(k);
                acgs(:, g_idx(k)) = double(ccg_3d(:, lid, lid));
            end
        end
    end
end
