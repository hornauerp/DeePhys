function diagnosticTrainLabels(ctc)
% DIAGNOSTICTRAINLABELS  2x3 diagnostic figure: "Are the training labels clean and separable?"
%
% Tiles:
%   (1,1) Filter funnel  — horizontal bar chart showing label-pipeline attrition
%   (1,2) PCA scatter    — training subset colored by pipeline stage
%   (1,3) Silhouette     — per-unit silhouette values for training set
%   (2,1) Cohen's d      — top 6 discriminative features
%   (2,2) CE distances   — counterexample distance distribution from inh centroid
%   (2,3) UMAP scatter   — unsupervised embedding colored by training assignment
%
% Each tile is guarded by try/catch.
% Calls ctc.saveDiagnosticFigure at the end.

arguments
    ctc CellTypeClassifier
end

C_INH   = [0.8 0.1 0.1];   % inhibitory: red
C_EXC   = [0.1 0.3 0.8];   % excitatory: blue
C_GREY  = [0.6 0.6 0.6];   % unlabeled:  grey
C_ORNG  = [0.9 0.5 0.0];   % geom-rejected: orange

fig = figure('Visible', 'on');
set(fig, 'Position', [100 100 1400 800]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Training labels diagnostics', 'FontWeight', 'bold');

% Gather shared objects (may be empty if called prematurely)
tl_struct = ctc.TrainLabels;
np        = ctc.NormalizationParams;
nf        = ctc.NormalizedFeatures;

% ── Helper: reconstruct normalized feature matrix for all unique units ────────
%   Mirrors the pipeline in classifyUnits exactly.
    function X_feat = buildXFeat()
        X = normalize(nf.X_pergroup, 'center', np.mu_global, 'scale', np.sigma_global);
        X(:, np.nan_cols) = [];
        X = X ./ np.scale;
        if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
            X = X(:, np.feature_selection_mask);
        end
        X_feat = X;
    end

% ── Tile (1,1): Filter funnel ─────────────────────────────────────────────────
nexttile(1);
try
    fc = tl_struct.filter_counts;
    inh_vals   = [fc.n_responsive_in, fc.n_after_outlier, fc.n_after_geometric];
    inh_labels = {'Candidates in', 'After outlier detection', 'After geom. filter'};
    exc_vals   = [fc.n_ce_in, fc.n_ce_after_outlier];
    exc_labels = {'CE candidates', 'After outlier detection'};

    n_inh = numel(inh_vals);
    n_exc = numel(exc_vals);

    % y positions: inh group at rows 1..n_inh, exc group offset by n_inh+1
    gap    = 1.5;
    y_inh  = (n_inh:-1:1)';
    y_exc  = (n_exc:-1:1)' + n_inh + gap;

    hold on;
    % Inhibitory bars
    for bi = 1:n_inh
        barh(y_inh(bi), inh_vals(bi), 'FaceColor', C_INH, 'EdgeColor', 'none');
        text(inh_vals(bi) + max(inh_vals)*0.01, y_inh(bi), ...
            num2str(inh_vals(bi)), 'VerticalAlignment', 'middle', 'FontSize', 8);
    end
    % Excitatory bars
    for bi = 1:n_exc
        barh(y_exc(bi), exc_vals(bi), 'FaceColor', C_EXC, 'EdgeColor', 'none');
        text(exc_vals(bi) + max([inh_vals, exc_vals])*0.01, y_exc(bi), ...
            num2str(exc_vals(bi)), 'VerticalAlignment', 'middle', 'FontSize', 8);
    end
    hold off;

    all_y      = [y_inh; y_exc];
    all_labels = [fliplr(inh_labels), fliplr(exc_labels)];
    yticks(all_y);
    yticklabels(all_labels);
    xlabel('Unit count');
    title('Label pipeline funnel');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (1,2): PCA scatter colored by pipeline stage ─────────────────────────
nexttile(2);
try
    X_feat = buildXFeat();

    % Subset rows
    X_sub  = X_feat(nf.subset_mask, :);
    n_comp = min(2, size(X_sub, 1) - 1);
    if n_comp < 2
        error('Too few training-subset units for PCA');
    end
    [~, scores] = pca(X_sub, 'NumComponents', 2);
    pc1 = scores(:, 1);
    pc2 = scores(:, 2);

    % Global unique indices for training units
    train_inh_global = tl_struct.sorted_train_ids(tl_struct.sorted_y_train == 2);
    train_exc_global = tl_struct.sorted_train_ids(tl_struct.sorted_y_train == 1);

    % Responsive candidates — local to subset (UnitID-based, order-independent)
    unit_ids_table   = string(ctc.FeatureStore.UnitTable.UnitID);
    resp_uids_diag   = unique(unit_ids_table(ctc.ResponsiveUnitIdx));
    unique_ud_ids_d  = string({nf.unique_ud.UnitID});  % row vector
    resp_unique_mask = ismember(unique_ud_ids_d, resp_uids_diag);
    subset_global    = find(nf.subset_mask);
    resp_local       = find(resp_unique_mask(nf.subset_mask));

    % Rejected-by-outlier local indices (local to subset)
    outlier_rejected_local = resp_local(tl_struct.resp_outlier_mask);

    % Rejected-by-geom local indices (local to subset, relative to after-outlier set)
    after_outlier_local = resp_local(~tl_struct.resp_outlier_mask);
    if ~isempty(tl_struct.resp_geom_mask) && ~isempty(after_outlier_local)
        geom_rejected_local = after_outlier_local(tl_struct.resp_geom_mask);
    else
        geom_rejected_local = zeros(1, 0);
    end

    % Map global unique indices to local-subset indices
    % subset_global(local) = global_unique_idx
    [~, inh_local] = ismember(train_inh_global, subset_global);
    [~, exc_local] = ismember(train_exc_global, subset_global);
    inh_local = inh_local(inh_local > 0);
    exc_local = exc_local(exc_local > 0);

    hold on;
    % All subset units: grey background
    scatter(pc1, pc2, 15, C_GREY, 'filled', 'MarkerFaceAlpha', 0.3, ...
        'DisplayName', 'Subset (unlabeled)');
    % Outlier-rejected inhibitory: black x
    if ~isempty(outlier_rejected_local)
        scatter(pc1(outlier_rejected_local), pc2(outlier_rejected_local), 30, ...
            'kx', 'LineWidth', 1.2, 'DisplayName', 'Inh: outlier rejected');
    end
    % Geom-rejected inhibitory: orange x
    if ~isempty(geom_rejected_local)
        scatter(pc1(geom_rejected_local), pc2(geom_rejected_local), 30, ...
            'x', 'MarkerEdgeColor', C_ORNG, 'LineWidth', 1.2, ...
            'DisplayName', 'Inh: geom rejected');
    end
    % Final inhibitory training units: red filled
    if ~isempty(inh_local)
        scatter(pc1(inh_local), pc2(inh_local), 40, C_INH, 'filled', ...
            'DisplayName', 'Inhibitory (train)');
    end
    % Final excitatory counterexamples: blue filled
    if ~isempty(exc_local)
        scatter(pc1(exc_local), pc2(exc_local), 40, C_EXC, 'filled', ...
            'DisplayName', 'Excitatory (train)');
    end
    hold off;
    xlabel('PC1');
    ylabel('PC2');
    title('Label pipeline (PC1 vs PC2)');
    legend('Location', 'best', 'Box', 'off', 'FontSize', 7);
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (1,3): Training-set silhouette ───────────────────────────────────────
nexttile(3);
try
    X_feat    = buildXFeat();
    train_idx = tl_struct.sorted_train_ids;   % global unique indices
    y_train   = tl_struct.sorted_y_train;

    X_train = X_feat(train_idx, :);
    n_train = size(X_train, 1);

    if n_train < 4
        error('Too few training units for silhouette');
    end

    % Reduce to at most 50 PCA components for speed
    n_pca = min([50, size(X_train, 1) - 1, size(X_train, 2)]);
    if n_pca >= 2
        [~, X_for_sil] = pca(X_train, 'NumComponents', n_pca);
    else
        X_for_sil = X_train;
    end

    sil = silhouette(X_for_sil, y_train');

    % Sort within each class and plot horizontal bars
    hold on;
    class_vals = unique(y_train);
    y_pos      = 0;
    for ci = 1:numel(class_vals)
        cv   = class_vals(ci);
        mask = y_train == cv;
        sv   = sort(sil(mask), 'ascend');
        n_c  = numel(sv);
        c    = [C_EXC; C_INH];
        col  = c(min(cv, 2), :);
        for ui = 1:n_c
            y_pos = y_pos + 1;
            barh(y_pos, sv(ui), 'FaceColor', col, 'EdgeColor', 'none');
        end
    end
    hold off;
    xline(0, '--k', 'LineWidth', 1);
    mean_sil = mean(sil);
    text(0.05, 0.95, sprintf('Mean = %.2f', mean_sil), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9);
    xlabel('Silhouette value');
    title(sprintf('Train silhouette (mean=%.2f)', mean_sil));
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (2,1): Top discriminative features (Cohen's d) ───────────────────────
nexttile(4);
try
    X_feat  = buildXFeat();
    train_idx = tl_struct.sorted_train_ids;
    y_train   = tl_struct.sorted_y_train;
    X_train   = X_feat(train_idx, :);

    mask1 = y_train == 1;
    mask2 = y_train == 2;
    mu1   = mean(X_train(mask1, :), 1);
    mu2   = mean(X_train(mask2, :), 1);
    s1    = std(X_train(mask1, :), 0, 1);
    s2    = std(X_train(mask2, :), 0, 1);
    n1    = sum(mask1);
    n2    = sum(mask2);

    pooled = sqrt(((n1 - 1) .* s1.^2 + (n2 - 1) .* s2.^2) ./ (n1 + n2 - 2));
    pooled(pooled == 0) = 1;
    d_vals = (mu2 - mu1) ./ pooled;   % positive = higher in inhibitory

    [~, sort_ord] = sort(abs(d_vals), 'descend');
    top6          = sort_ord(1:min(6, numel(sort_ord)));
    d_top         = d_vals(top6);

    feat_names = np.feat_names_trimmed;
    fn_top     = feat_names(top6);
    % Shorten names to last 20 chars
    for fi = 1:numel(fn_top)
        s = char(fn_top(fi));
        if numel(s) > 20
            fn_top(fi) = string(s(end-19:end));
        end
    end

    hold on;
    for fi = 1:numel(d_top)
        col = C_INH;
        if d_top(fi) < 0
            col = C_EXC;
        end
        barh(fi, d_top(fi), 'FaceColor', col, 'EdgeColor', 'none');
    end
    hold off;
    xline(0, '--k', 'LineWidth', 1);
    yticks(1:numel(d_top));
    yticklabels(fn_top);
    xlabel("Cohen's d (inh - exc)");
    title('Top 6 discriminative features');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (2,2): Counterexample distance distribution ─────────────────────────
nexttile(5);
try
    dists = tl_struct.ce_distances_from_inh_centroid;
    if isempty(dists)
        error('No distance data available');
    end

    dists = dists(:);

    % Inhibitory training unit distances from inh centroid
    X_feat       = buildXFeat();
    train_idx    = tl_struct.sorted_train_ids;
    y_train      = tl_struct.sorted_y_train;
    inh_centroid = tl_struct.inh_centroid;
    X_inh_train  = X_feat(train_idx(y_train == 2), :);

    if ~isempty(X_inh_train) && ~isempty(inh_centroid)
        inh_dists = pdist2(inh_centroid, X_inh_train, 'correlation');
        inh_dists = inh_dists(:);
    else
        inh_dists = [];
    end

    hold on;
    histogram(dists, 20, 'FaceColor', [0.4 0.7 1.0], 'FaceAlpha', 0.6, ...
        'EdgeColor', 'none', 'DisplayName', 'CE pool');

    if ~isempty(inh_dists)
        histogram(inh_dists, 20, 'FaceColor', C_INH, 'FaceAlpha', 0.5, ...
            'EdgeColor', 'none', 'DisplayName', 'Inh training');
    end

    % Vertical line at selection threshold (n_pick-th value)
    fc     = tl_struct.filter_counts;
    n_pick = fc.n_ce_in;   % number of CE candidates initially selected
    if n_pick > 0 && n_pick <= numel(dists)
        sorted_d  = sort(dists, 'descend');
        thresh_d  = sorted_d(n_pick);
        xline(thresh_d, '--r', 'LineWidth', 1.5, 'DisplayName', 'Selection threshold');
    end
    hold off;
    xlabel('Correlation distance from inh centroid');
    ylabel('Count');
    title('CE distance distribution');
    legend('Location', 'best', 'Box', 'off', 'FontSize', 8);
    box off;
catch ME
    cla; axis off;
    if contains(ME.message, 'No distance')
        text(0.5, 0.5, 'No distance data available', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 9);
    else
        text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', 8, 'Color', [0.6 0 0]);
    end
end

% ── Tile (2,3): Unsupervised UMAP colored by training assignment ───────────────
nexttile(6);
try
    unsup = ctc.Reduction.Unsupervised;
    if isempty(unsup)
        error('UMAP not available');
    end

    % unsup rows correspond to training-subset unique units (subset_mask == true)
    n_sub = size(unsup, 1);
    subset_global = find(nf.subset_mask);

    % Determine per-subset-row assignment
    train_inh_global = tl_struct.sorted_train_ids(tl_struct.sorted_y_train == 2);
    train_exc_global = tl_struct.sorted_train_ids(tl_struct.sorted_y_train == 1);

    [~, inh_local] = ismember(train_inh_global, subset_global);
    [~, exc_local] = ismember(train_exc_global, subset_global);
    inh_local = inh_local(inh_local > 0);
    exc_local = exc_local(exc_local > 0);

    unlabeled_mask = true(n_sub, 1);
    unlabeled_mask(inh_local) = false;
    unlabeled_mask(exc_local) = false;

    hold on;
    if any(unlabeled_mask)
        scatter(unsup(unlabeled_mask, 1), unsup(unlabeled_mask, 2), ...
            10, C_GREY, 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Unlabeled');
    end
    if ~isempty(exc_local)
        scatter(unsup(exc_local, 1), unsup(exc_local, 2), ...
            40, C_EXC, 'filled', 'DisplayName', 'Excitatory (train)');
    end
    if ~isempty(inh_local)
        scatter(unsup(inh_local, 1), unsup(inh_local, 2), ...
            40, C_INH, 'filled', 'DisplayName', 'Inhibitory (train)');
    end
    hold off;
    xlabel('UMAP 1');
    ylabel('UMAP 2');
    title('Unsupervised UMAP (training labels)');
    legend('Location', 'best', 'Box', 'off', 'FontSize', 8);
    box off;
catch ME
    cla; axis off;
    if contains(ME.message, 'UMAP not available')
        text(0.5, 0.5, 'UMAP not available', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 9);
    else
        text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', 8, 'Color', [0.6 0 0]);
    end
end

% ── Save ──────────────────────────────────────────────────────────────────────
ctc.saveDiagnosticFigure(fig, 'train_labels');

end
