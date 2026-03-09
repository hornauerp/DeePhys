function plotCellTypeFeatures(unit_list, unit_labels, train_labels)
% PLOTCELLTYPEFEATURES  Waveform and ACG heatmaps split by cell type and train/test.
%
% Produces an 8-panel (2 rows x 4 columns) figure:
%   Rows:    Excitatory (top) | Inhibitory (bottom)
%   Columns: Train waveforms | Train ACGs | Test waveforms | Test ACGs
%
% Waveforms are sorted by trough position; ACGs by peak position.
%
% INPUTS:
%   unit_list    - (1xN) Unit array (ctc.UnitList)
%   unit_labels  - (1xN) double: 1=excitatory, 2=inhibitory, NaN=unclassified
%   train_labels - ctc.TrainLabels struct with fields umap_train_idx / umap_test_idx

train_mask = train_labels.umap_train_idx(:)';
test_mask  = train_labels.umap_test_idx(:)';

% Feature matrices: units x samples/lags, row-normalised
all_wf   = horzcat(unit_list.ReferenceWaveform)';
all_acgs = horzcat(unit_list.FullACG)';

wf_scale  = max(abs(all_wf),  [], 2); wf_scale(wf_scale   == 0) = 1;
acg_scale = max(all_acgs,     [], 2); acg_scale(acg_scale  == 0) = 1;
norm_wf   = all_wf   ./ wf_scale;
norm_acgs = all_acgs ./ acg_scale;

% Layout: 2 rows (exc/inh) x 4 cols (train-wf | train-acg | test-wf | test-acg)
cell_type_vals = [1, 2];
type_labels    = ["Excitatory", "Inhibitory"];
split_masks    = {train_mask, train_mask, test_mask, test_mask};
col_labels     = ["Train — waveforms", "Train — ACGs", "Test — waveforms", "Test — ACGs"];
use_acg        = [false, true, false, true];

figure('Color', 'w');
tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Cell-type classification — waveforms & ACGs (top: excitatory | bottom: inhibitory)', ...
    'FontWeight', 'normal');

for r = 1:2
    ct_mask = unit_labels == cell_type_vals(r);

    for c = 1:4
        m = ct_mask & split_masks{c};

        if use_acg(c)
            sub = norm_acgs(m, :);
            if ~isempty(sub)
                [~, peak_pos]  = max(sub, [], 2);
                [~, ord]       = sort(peak_pos);
                sub            = sub(ord, :);
            end
        else
            sub = norm_wf(m, :);
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
            imagesc(sub); colorbar; axis tight
            ylabel(sprintf('n = %i', size(sub, 1)));
            if use_acg(c); xlabel('Lag (bins)'); else; xlabel('Sample'); end
        end

        ttl = col_labels(c);
        if c == 1
            ttl = type_labels(r) + "  |  " + ttl;
        end
        title(ttl, 'FontWeight', 'normal');
    end
end
