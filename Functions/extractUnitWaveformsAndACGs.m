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
unit_srs = arrayfun(@(u) u.MEArecording.RecordingInfo.SamplingRate, unit_list);
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
    use_cached = (n_bins_cached == n_bins_harm);
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
    fprintf('Harmonization ACG params (%i bins) differ from cached %s (%i bins) — recomputing per recording\n', ...
        n_bins_harm, acg_source, n_bins_cached);

    % Group units by parent MEArecording for one batch CCG call per recording.
    % All spike times are already in memory — no .npy files are re-read.
    rec_list   = {};
    rec_groups = {};
    for u = 1:N
        rec   = unit_list(u).MEArecording;
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
        rec   = rec_list{r};
        g_idx = rec_groups{r};          % global indices into unit_list
        sr_r  = rec.RecordingInfo.SamplingRate;

        % Build concatenated spike train; track compact local unit index
        all_times    = [];
        all_unit_ids = [];
        local_ids    = zeros(1, numel(g_idx));  % local CCG index per unit (0 = no spikes)
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
            continue  % all units in this recording have no spikes
        end

        [~, sort_idx]  = sort(all_times);
        all_times      = all_times(sort_idx);
        all_unit_ids   = all_unit_ids(sort_idx);

        [ccg_3d, ~] = CCG(all_times, all_unit_ids, ...
            'binSize', ph.ACGBinSize, ...
            'duration', 2 * ph.ACGLag, ...
            'Fs', 1 / sr_r);

        % Extract ACG (diagonal) for each unit in this recording
        for k = 1:numel(g_idx)
            if local_ids(k) > 0
                lid = local_ids(k);
                acgs(:, g_idx(k)) = double(ccg_3d(:, lid, lid));
            end
        end
    end
end
end
