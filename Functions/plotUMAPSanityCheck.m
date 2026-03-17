function plotUMAPSanityCheck(ctc)
% PLOTUMAPSANITYCHECK  Scatter-plot UMAP embeddings to sanity-check classification.
%
% Produces a figure with 2 or 3 tiles depending on what is available:
%
%   Tile 1 (bootstrap path only):
%     Unsupervised UMAP — all units in the training subset.
%     Shows: non-responsive (grey), responsive outliers removed by iforest
%     (black x), clean inhibitory candidates (red), excitatory
%     counterexamples (blue).
%
%   Tile 2 (bootstrap path only):
%     Unsupervised UMAP — same embedding, coloured by final predicted
%     cell-type labels (excitatory/inhibitory/unclassified).
%
%   Tile 3:
%     Supervised UMAP — training units coloured by known label.
%
%   Tile 4:
%     Supervised UMAP — test units coloured by predicted label.
%     Unclassified units (NaN label) are shown in grey.
%
% Call after ctc.classifyUnits().
% Requires ctc.Reduction.Train and ctc.Reduction.Test to be set.
% Tile 1 additionally requires ctc.Reduction.Unsupervised and
% ctc.ResponsiveUnitIdx (set by identifyResponsiveUnits / generateTrainLabels).

assert(~isempty(ctc.Reduction.Train), ...
    'ctc.Reduction.Train is empty — run classifyUnits() first.');

exc_color = [0.22, 0.49, 0.87];  % blue  — excitatory
inh_color = [0.87, 0.28, 0.20];  % red   — inhibitory
gray      = [0.70, 0.70, 0.70];
black     = [0.15, 0.15, 0.15];
colors    = [exc_color; inh_color];

has_unsup = ~isempty(ctc.Reduction.Unsupervised) && ~isempty(ctc.ResponsiveUnitIdx);
n_tiles   = 2 + 2 * has_unsup;  % 2 extra tiles when unsupervised data exists

figure('Color', 'w');
tl = tiledlayout(1, n_tiles, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'UMAP sanity check', 'FontWeight', 'normal');

% --- Tile 1: unsupervised embedding with outlier/candidate detail ---
if has_unsup
    r    = ctc.Reduction.Unsupervised;
    resp = ctc.ResponsiveUnitIdx(:)';
    N    = numel(ctc.UnitList);

    % Determine which subset was used for UMAP (matches generateTrainLabels logic)
    p_umap = ctc.Parameters.UMAP;
    if ~isempty(p_umap.TrainingCultureIdx)
        subset_mask = buildCultureUnitMask(ctc, p_umap.TrainingCultureIdx);
    else
        subset_mask = true(1, N);
    end

    % Build per-unit category within the subset
    % Categories: 1=non-responsive, 2=outlier, 3=inhibitory candidate, 4=excitatory candidate
    subset_global = find(subset_mask);
    cat = ones(1, size(r, 1));  % default: non-responsive

    % Responsive candidates (before outlier removal)
    subset_responsive = resp(subset_mask);
    cat(subset_responsive) = 3;  % inhibitory candidate (clean)

    % Outliers: responsive candidates removed by isolation forest
    if isfield(ctc.TrainLabels, 'outlier_global_idx')
        outlier_global = ctc.TrainLabels.outlier_global_idx;
        [~, outlier_local] = ismember(outlier_global, subset_global);
        outlier_local(outlier_local == 0) = [];
        cat(outlier_local) = 2;  % outlier
    end

    % Excitatory counterexamples
    if isfield(ctc.TrainLabels, 'excitatory_candidate_global_idx')
        exc_global = ctc.TrainLabels.excitatory_candidate_global_idx;
        [~, exc_local] = ismember(exc_global, subset_global);
        exc_local(exc_local == 0) = [];
        cat(exc_local) = 4;  % excitatory candidate
    end

    nexttile; hold on;

    % Layer 1: non-responsive (grey, background)
    m1 = cat == 1;
    if any(m1)
        scatter(r(m1, 1), r(m1, 2), 6, gray, 'filled', ...
            'MarkerFaceAlpha', 0.2, 'DisplayName', sprintf('Non-responsive (n=%d)', sum(m1)));
    end

    % Layer 2: excitatory counterexamples (blue)
    m4 = cat == 4;
    if any(m4)
        scatter(r(m4, 1), r(m4, 2), 12, exc_color, 'filled', ...
            'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('Excitatory candidate (n=%d)', sum(m4)));
    end

    % Layer 3: clean inhibitory candidates (red)
    m3 = cat == 3;
    if any(m3)
        scatter(r(m3, 1), r(m3, 2), 12, inh_color, 'filled', ...
            'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('Inhibitory candidate (n=%d)', sum(m3)));
    end

    % Layer 4: outliers (black x, on top)
    m2 = cat == 2;
    if any(m2)
        scatter(r(m2, 1), r(m2, 2), 18, black, 'x', ...
            'LineWidth', 1, 'DisplayName', sprintf('Outlier removed (n=%d)', sum(m2)));
    end

    xlabel('UMAP 1'); ylabel('UMAP 2');
    title('Unsupervised — label generation', 'FontWeight', 'normal');
    legend('Location', 'best', 'FontSize', 7, 'Box', 'off');
    box off; hold off;

    % --- Tile 2: unsupervised embedding coloured by predicted cell type ---
    nexttile; hold on;

    % Map global UnitLabels onto the subset used for unsupervised UMAP
    subset_labels = ctc.UnitLabels(subset_mask);

    m_exc = subset_labels == 1;
    m_inh = subset_labels == 2;
    m_unc = isnan(subset_labels) | (subset_labels ~= 1 & subset_labels ~= 2);

    if any(m_unc)
        scatter(r(m_unc, 1), r(m_unc, 2), 6, gray, 'filled', ...
            'MarkerFaceAlpha', 0.2, 'DisplayName', sprintf('Unclassified (n=%d)', sum(m_unc)));
    end
    if any(m_exc)
        scatter(r(m_exc, 1), r(m_exc, 2), 10, exc_color, 'filled', ...
            'MarkerFaceAlpha', 0.5, 'DisplayName', sprintf('Excitatory (n=%d)', sum(m_exc)));
    end
    if any(m_inh)
        scatter(r(m_inh, 1), r(m_inh, 2), 10, inh_color, 'filled', ...
            'MarkerFaceAlpha', 0.5, 'DisplayName', sprintf('Inhibitory (n=%d)', sum(m_inh)));
    end

    xlabel('UMAP 1'); ylabel('UMAP 2');
    title('Unsupervised — predicted labels', 'FontWeight', 'normal');
    legend('Location', 'best', 'FontSize', 7, 'Box', 'off');
    box off; hold off;
end

% --- Tile 3: supervised embedding — training units ---
r_train  = ctc.Reduction.Train;
y_train  = ctc.TrainLabels.sorted_y_train(:)';

nexttile; hold on;
for cl = [1, 2]
    m = y_train == cl;
    scatter(r_train(m, 1), r_train(m, 2), 10, colors(cl, :), 'filled', ...
        'MarkerFaceAlpha', 0.7, 'DisplayName', labelName(cl));
end
xlabel('UMAP 1'); ylabel('UMAP 2');
title('Supervised — training set', 'FontWeight', 'normal');
legend('Location', 'best', 'FontSize', 7);
box off; hold off;

% --- Tile 4: supervised embedding — test units ---
r_test       = ctc.Reduction.Test;
test_mask    = ctc.TrainLabels.umap_test_idx(:)';
labels_test  = ctc.UnitLabels(test_mask);

nexttile; hold on;
for cl = [1, 2]
    m = labels_test == cl;
    scatter(r_test(m, 1), r_test(m, 2), 6, colors(cl, :), 'filled', ...
        'MarkerFaceAlpha', 0.4, 'DisplayName', labelName(cl));
end
m_nan = isnan(labels_test);
if any(m_nan)
    scatter(r_test(m_nan, 1), r_test(m_nan, 2), 6, gray, 'filled', ...
        'MarkerFaceAlpha', 0.3, 'DisplayName', 'Unclassified');
end
xlabel('UMAP 1'); ylabel('UMAP 2');
title('Supervised — test set (predicted)', 'FontWeight', 'normal');
legend('Location', 'best', 'FontSize', 7);
box off; hold off;

end

function name = labelName(cl)
if cl == 1; name = 'Excitatory'; else; name = 'Inhibitory'; end
end

function mask = buildCultureUnitMask(ctc, culture_indices)
% Build a logical mask over ctc.UnitList selecting units from specified cultures.
    rg = ctc.RecordingGroup;
    mask = false(1, numel(ctc.UnitList));
    offset = 0;
    for c = 1:length(rg.Cultures)
        culture = rg.Cultures(c);
        n_units_c = numel(culture.Units);
        if ismember(c, culture_indices)
            mask(offset+1 : offset+n_units_c) = true;
        end
        offset = offset + n_units_c;
    end
end
