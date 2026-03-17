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
train_mask = ctc.TrainLabels.umap_train_idx(:)';
test_mask  = ctc.TrainLabels.umap_test_idx(:)';
labels     = ctc.UnitLabels;

% Harmonized data: (N_samples × N_units) — transpose to (N_units × N_samples)
all_wf   = ctc.HarmonizedWaveforms';
all_acgs = ctc.HarmonizedACGs';

% Row-normalise for heatmap display
wf_scale  = max(abs(all_wf),  [], 2); wf_scale(wf_scale   == 0) = 1;
acg_scale = max(all_acgs,     [], 2); acg_scale(acg_scale  == 0) = 1;
norm_wf   = all_wf   ./ wf_scale;
norm_acgs = all_acgs ./ acg_scale;

% Time axes for x-labels
ph    = ctc.Parameters.Harmonization;
t_wf  = (0:size(all_wf, 2)-1) / ph.WaveformTargetSamplingRate * 1000;  % ms
t_acg = linspace(-ph.ACGLag, ph.ACGLag, size(all_acgs, 2)) * 1000;    % ms

% ── Figure 1: heatmap grid ─────────────────────────────────────────────────
cell_type_vals = [1, 2];
type_labels    = ["Excitatory", "Inhibitory"];
split_masks    = {train_mask, train_mask, test_mask, test_mask};
col_labels     = ["Train - waveforms", "Train - ACGs", "Test - waveforms", "Test - ACGs"];
use_acg        = [false, true, false, true];

figure('Color', 'w', 'Name', 'Cell-type features (harmonized)');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Harmonized classifier input (top: excitatory | bottom: inhibitory)', ...
    'FontWeight', 'normal');

for r = 1:2
    ct_mask = labels == cell_type_vals(r);

    for c = 1:4
        m = ct_mask & split_masks{c};

        if use_acg(c)
            sub = norm_acgs(m, :);
            t_x = t_acg;
            x_label = 'Lag (ms)';
            if ~isempty(sub)
                [~, peak_pos] = max(sub, [], 2);
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
colors = struct('exc', [0.2 0.4 0.8], 'inh', [0.8 0.2 0.2]);

figure('Color', 'w', 'Name', 'Average waveforms by cell type');
tl2 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Mean waveform (+/- SEM) by cell type', 'FontWeight', 'normal');

split_names  = ["Train", "Test"];
split_idx    = {train_mask, test_mask};

for s = 1:2
    nexttile; hold on;
    mask_s = split_idx{s};

    for ct = 1:2
        m = mask_s & (labels == cell_type_vals(ct));
        wf_sub = all_wf(m, :);  % unnormalised — preserves amplitude
        if isempty(wf_sub); continue; end

        mu  = mean(wf_sub, 1);
        sem = std(wf_sub, 0, 1) / sqrt(size(wf_sub, 1));

        if ct == 1
            clr = colors.exc;
        else
            clr = colors.inh;
        end

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
