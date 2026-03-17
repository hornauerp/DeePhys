function plotCellTypeFeatures(ctc)
% PLOTCELLTYPEFEATURES  Waveform and ACG heatmaps split by cell type and train/test.
%
% Uses the harmonized data stored on ctc (HarmonizedWaveforms, HarmonizedACGs)
% so the plots reflect the actual classifier input — upsampled waveforms and
% ACGs at the configured Harmonization parameters.
%
% Produces two figures:
%   Figure 1 — 2×4 heatmap grid:
%     Rows:    Excitatory (top) | Inhibitory (bottom)
%     Columns: Train waveforms | Train ACGs | Test waveforms | Test ACGs
%     Waveforms sorted by trough position; ACGs sorted by peak position.
%
%   Figure 2 — 2×2 average waveform overlay:
%     Rows:    Train (top) | Test (bottom)
%     Columns: Overlay with mean+SEM (left) | Normalised overlay (right)
%
% INPUT:
%   ctc - CellTypeClassifier with HarmonizedWaveforms, HarmonizedACGs,
%         TrainLabels, and UnitLabels populated (run classifyUnits first)

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.HarmonizedWaveforms), ...
    'ctc.HarmonizedWaveforms is empty — run classifyUnits() or classifyUnitsWithExternalTrain() first.');
assert(~isempty(ctc.TrainLabels), ...
    'ctc.TrainLabels is empty — run generateTrainLabels() first.');

% ── Unpack ──────────────────────────────────────────────────────────────────
labels = ctc.UnitLabels;
N      = numel(labels);
ph     = ctc.Parameters.Harmonization;

% Detect transfer mode: TrainLabels indices refer to a different (source)
% unit list, so train_mask/test_mask don't apply to this ctc's data.
% In transfer mode, all units in this ctc are "test" (target) units.
is_transfer = false;
if isfield(ctc.TrainLabels, 'umap_train_idx')
    train_idx_raw = ctc.TrainLabels.umap_train_idx(:)';
    test_idx_raw  = ctc.TrainLabels.umap_test_idx(:)';
    % Check if the masks are compatible with this ctc's unit count
    if max([find(train_idx_raw, 1, 'last'), find(test_idx_raw, 1, 'last')]) <= N ...
            && numel(train_idx_raw) == N
        train_mask = logical(train_idx_raw);
        test_mask  = logical(test_idx_raw);
    else
        is_transfer = true;
    end
else
    is_transfer = true;
end

if is_transfer
    train_mask = false(1, N);      % no training units in target ctc
    test_mask  = true(1, N);       % all units are test
end

% Harmonized data: (N_samples × N_units) — transpose to (N_units × N_samples)
all_wf   = ctc.HarmonizedWaveforms';
all_acgs = ctc.HarmonizedACGs';

% Row-normalise for heatmap display
wf_scale  = max(abs(all_wf),  [], 2); wf_scale(wf_scale   == 0) = 1;
acg_scale = max(all_acgs,     [], 2); acg_scale(acg_scale  == 0) = 1;
norm_wf   = all_wf   ./ wf_scale;
norm_acgs = all_acgs ./ acg_scale;

% ── Time axes ──────────────────────────────────────────────────────────────
% Waveform axis: trough-centred using Harmonization time window
n_wf_samp = size(all_wf, 2);
if isfield(ph, 'WaveformPreTrough') && isfield(ph, 'WaveformPostTrough')
    t_wf = linspace(-ph.WaveformPreTrough, ph.WaveformPostTrough, n_wf_samp);  % ms
else
    % Legacy fallback: assume zero-based
    wf_sr = ph.WaveformTargetSamplingRate;
    t_wf  = (0:n_wf_samp-1) / wf_sr * 1000;
end

% ACG axis: symmetric around zero
n_acg_bins = size(all_acgs, 2);
half_bins  = floor(n_acg_bins / 2);
t_acg = (-half_bins:half_bins) * ph.ACGBinSize * 1000;  % ms
if numel(t_acg) ~= n_acg_bins
    t_acg = linspace(-ph.ACGLag, ph.ACGLag, n_acg_bins) * 1000;
end

% ── Figure 1: heatmap grid ─────────────────────────────────────────────────
cell_type_vals = [1, 2];
type_labels    = ["Excitatory", "Inhibitory"];

if is_transfer
    % Transfer mode: 2×2 grid (target units only, no train/test split)
    split_masks = {test_mask, test_mask};
    col_labels  = ["Waveforms", "ACGs"];
    use_acg     = [false, true];
    n_cols      = 2;
    fig_title   = 'Transfer classifier — target units';
else
    % Standard mode: 2×4 grid (train + test)
    split_masks = {train_mask, train_mask, test_mask, test_mask};
    col_labels  = ["Train - waveforms", "Train - ACGs", "Test - waveforms", "Test - ACGs"];
    use_acg     = [false, true, false, true];
    n_cols      = 4;
    fig_title   = 'Harmonized classifier input (top: excitatory | bottom: inhibitory)';
end

figure('Color', 'w', 'Name', 'Cell-type features (harmonized)');
tl = tiledlayout(2, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, fig_title, 'FontWeight', 'normal');

for r = 1:2
    ct_mask = labels == cell_type_vals(r);

    for c = 1:n_cols
        m = ct_mask & split_masks{c};

        if use_acg(c)
            sub = norm_acgs(m, :);
            t_x = t_acg;
            x_label = 'Lag (ms)';
            if ~isempty(sub)
                half = floor(size(sub, 2) / 2);
                [~, peak_pos] = max(sub(:, 1:half), [], 2);
                [~, ord]      = sort(peak_pos);
                sub           = sub(ord, :);
            end
        else
            sub = norm_wf(m, :);
            t_x = t_wf;
            x_label = 'Time (ms)';
            if ~isempty(sub)
                [~, trough_pos] = min(sub, [], 2);
                [~, ord]        = sort(trough_pos);
                sub             = sub(ord, :);
            end
        end

        nexttile;
        if isempty(sub)
            text(0.5, 0.5, 'no units', 'HorizontalAlignment', 'center', ...
                'Units', 'normalized');
            axis off
        else
            imagesc(t_x, 1:size(sub,1), sub); colorbar; axis tight
            if ~use_acg(c)
                clim([-1 0.5]);
            end
            xlabel(x_label);
            ylabel(sprintf('n = %i', size(sub, 1)));
        end

        ttl = col_labels(c);
        if c == 1
            ttl = type_labels(r) + "  |  " + ttl;
        end
        title(ttl, 'FontWeight', 'normal');
    end
end

% ── Figure 2: average waveforms per cell type ──────────────────────────────
clrs = struct('exc', [0.2 0.4 0.8], 'inh', [0.8 0.2 0.2]);

if is_transfer
    split_names = "Target";
    split_idx   = {test_mask};
else
    split_names = ["Train", "Test"];
    split_idx   = {train_mask, test_mask};
end

figure('Color', 'w', 'Name', 'Average waveforms by cell type');
tl2 = tiledlayout(1, numel(split_idx), 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Mean waveform (\pm SEM) by cell type', 'FontWeight', 'normal');

for s = 1:numel(split_idx)
    nexttile; hold on;
    mask_s = split_idx{s};

    for ct = 1:2
        m = mask_s & (labels == cell_type_vals(ct));
        wf_sub = all_wf(m, :);
        if isempty(wf_sub); continue; end

        mu  = mean(wf_sub, 1);
        sem = std(wf_sub, 0, 1) / sqrt(size(wf_sub, 1));

        if ct == 1; clr = clrs.exc; else; clr = clrs.inh; end

        fill([t_wf, fliplr(t_wf)], [mu+sem, fliplr(mu-sem)], clr, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(t_wf, mu, 'Color', clr, 'LineWidth', 1.5);
    end

    ylim([-1 0.5]);
    xlabel('Time (ms)');
    ylabel('Amplitude (norm.)');
    title(split_names(s), 'FontWeight', 'normal');
    legend({'', 'Excitatory', '', 'Inhibitory'}, 'Location', 'best', 'Box', 'off');
    box off;
    hold off;
end
end
