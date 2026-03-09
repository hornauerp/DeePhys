function plotUMAPSanityCheck(ctc)
% PLOTUMAPSANITYCHECK  Scatter-plot UMAP embeddings to sanity-check classification.
%
% Produces a figure with 2 or 3 tiles depending on what is available:
%
%   Tile 1 (bootstrap path only):
%     Unsupervised UMAP — all units.
%     Responsive (inhibitory candidate) units are highlighted in red;
%     non-responsive in grey.
%
%   Tile 2:
%     Supervised UMAP — training units coloured by known label.
%
%   Tile 3:
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
colors    = [exc_color; inh_color];

has_unsup = ~isempty(ctc.Reduction.Unsupervised) && ~isempty(ctc.ResponsiveUnitIdx);
n_tiles   = 2 + has_unsup;

figure('Color', 'w');
tl = tiledlayout(1, n_tiles, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'UMAP sanity check', 'FontWeight', 'normal');

% --- Tile 1: unsupervised embedding (bootstrap path only) ---
if has_unsup
    r    = ctc.Reduction.Unsupervised;
    resp = ctc.ResponsiveUnitIdx(:)';

    nexttile; hold on;
    scatter(r(~resp, 1), r(~resp, 2), 6, gray, 'filled', ...
        'MarkerFaceAlpha', 0.3, 'DisplayName', 'Non-responsive');
    scatter(r(resp, 1),  r(resp, 2),  10, inh_color, 'filled', ...
        'MarkerFaceAlpha', 0.7, 'DisplayName', 'Responsive candidate');
    xlabel('UMAP 1'); ylabel('UMAP 2');
    title('Unsupervised — responsive candidates', 'FontWeight', 'normal');
    legend('Location', 'best', 'FontSize', 7);
    box off; hold off;
end

% --- Tile 2: supervised embedding — training units ---
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

% --- Tile 3: supervised embedding — test units ---
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
