function metrics = validateClassification(pred_labels, true_labels)
% VALIDATECLASSIFICATION  Confusion matrix and performance metrics vs ground truth.
%
% Evaluates cell-type classifier performance against an independent ground
% truth.  Units where true_labels is NaN (unknown ground truth) are silently
% excluded; units where pred_labels is NaN (unclassified) are counted
% separately but also excluded from metric computation.
%
% INPUTS:
%   pred_labels - (1xN) predicted labels: 1=excitatory, 2=inhibitory, NaN=unclassified
%   true_labels - (1xN) ground truth labels: 1=excitatory, 2=inhibitory, NaN=unknown
%
% OUTPUT:
%   metrics - struct with Accuracy, Precision, Recall, F1score, Specificity
%             per class and macro-averaged (from multiclass_metrics_common)

pred_labels = pred_labels(:)';
true_labels = true_labels(:)';

has_gt   = ~isnan(true_labels);
has_pred = ~isnan(pred_labels);
valid    = has_gt & has_pred;

n_gt      = sum(has_gt);
n_valid   = sum(valid);
n_unclass = sum(has_gt & ~has_pred);

fprintf('--- Validation summary ---\n');
fprintf('Units with ground truth:      %i\n', n_gt);
fprintf('  classified (evaluated):     %i  (%.1f%%)\n', n_valid,   100*n_valid/max(n_gt,1));
fprintf('  unclassified (excluded):    %i  (%.1f%%)\n', n_unclass, 100*n_unclass/max(n_gt,1));

y_true = true_labels(valid);
y_pred = pred_labels(valid);

fprintf('True label breakdown:  %i excitatory | %i inhibitory\n', ...
    sum(y_true == 1), sum(y_true == 2));
fprintf('Pred label breakdown:  %i excitatory | %i inhibitory\n', ...
    sum(y_pred == 1), sum(y_pred == 2));

% 2x2 confusion matrix: rows = true, cols = predicted, order [exc=1, inh=2]
C = confusionmat(y_true, y_pred, 'Order', [1, 2]);

% Metrics via existing helper (handles both 2-class and N-class)
metrics = multiclass_metrics_common(C);

fprintf('\nAccuracy:    %.3f\n', metrics.Accuracy);
fprintf('Precision:   %.3f  (macro)\n', metrics.Precision);
fprintf('Recall:      %.3f  (macro)\n', metrics.Recall);
fprintf('F1 score:    %.3f  (macro)\n', metrics.F1score);
fprintf('Specificity: %.3f\n',          metrics.Specificity);

plotConfusionMatrixPair(C);
end

% -------------------------------------------------------------------------
function plotConfusionMatrixPair(C)

class_names = {'Excitatory', 'Inhibitory'};
n           = size(C, 1);
row_sums    = sum(C, 2);
C_pct       = 100 * C ./ max(row_sums, 1);   % row-normalised, avoid /0

% White-to-blue colormap (base MATLAB, no extra toolbox needed)
blues = [linspace(1, 0.08, 64)', linspace(1, 0.38, 64)', ones(64, 1)];

figure('Color', 'w');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Classification validation — confusion matrix', 'FontWeight', 'normal');

% --- Tile 1: raw counts ---
nexttile;
imagesc(C);
colormap(gca, blues); colorbar;
axis square;
set(gca, 'XTick', 1:n, 'XTickLabel', class_names, ...
         'YTick', 1:n, 'YTickLabel', class_names, ...
         'TickLength', [0 0]);
xlabel('Predicted'); ylabel('True');
title('Count', 'FontWeight', 'normal');
annotateMatrix(C, @(v) num2str(v), max(C(:)));

% --- Tile 2: row-normalised (%) ---
nexttile;
imagesc(C_pct, [0 100]);
colormap(gca, blues); colorbar;
axis square;
set(gca, 'XTick', 1:n, 'XTickLabel', class_names, ...
         'YTick', 1:n, 'YTickLabel', class_names, ...
         'TickLength', [0 0]);
xlabel('Predicted'); ylabel('True');
title('Row-normalised (%)', 'FontWeight', 'normal');
annotateMatrix(C_pct, @(v) sprintf('%.1f%%', v), 100);
end

% -------------------------------------------------------------------------
function annotateMatrix(M, fmtFn, maxVal)
[nr, nc] = size(M);
for r = 1:nr
    for c = 1:nc
        col = pickTextColor(M(r,c), maxVal);
        text(c, r, fmtFn(M(r,c)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontWeight', 'bold', 'Color', col);
    end
end
end

% -------------------------------------------------------------------------
function col = pickTextColor(val, maxVal)
% White text on dark cells, black on light
if val > 0.55 * maxVal
    col = [1 1 1];
else
    col = [0 0 0];
end
end
