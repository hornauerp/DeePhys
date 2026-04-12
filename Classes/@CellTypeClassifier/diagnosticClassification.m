function diagnosticClassification(ctc)
% DIAGNOSTICCLASSIFICATION  2x3 diagnostic figure: "Is the classification biologically plausible and confident?"
%
% Tiles:
%   (1,1) Confidence distribution   — histogram of kNN confidence for test units
%   (1,2) Per-chip I/E fraction      — bar chart with expected range shading
%   (1,3) UMAP test confidence       — supervised embedding colored by confidence
%   (2,1) Mean waveform by class     — mean +/- SEM with peak-to-trough annotation
%   (2,2) Mean ACG by class          — mean +/- SEM
%   (2,3) Feature effect sizes       — Cohen's d from training vs prediction
%
% Each tile is guarded by try/catch.
% Calls ctc.saveDiagnosticFigure at the end.

arguments
    ctc CellTypeClassifier
end

C_INH  = [0.8 0.1 0.1];   % inhibitory: red
C_EXC  = [0.1 0.3 0.8];   % excitatory: blue

fig = figure('Visible', 'on');
set(fig, 'Position', [100 100 1400 800]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Classification diagnostics', 'FontWeight', 'bold');

% ── Helper: reconstruct normalized feature matrix for all unique units ────────
    function X_feat = buildXFeat()
        np = ctc.NormalizationParams;
        nf = ctc.NormalizedFeatures;
        X  = normalize(nf.X_pergroup, 'center', np.mu_global, 'scale', np.sigma_global);
        X(:, np.nan_cols) = [];
        X = X ./ np.scale;
        if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
            X = X(:, np.feature_selection_mask);
        end
        X_feat = X;
    end

% ── Helper: identify test-unit mask (full array, excluding training units) ────
    function test_mask = getTestMask()
        nf     = ctc.NormalizedFeatures;
        tl_s   = ctc.TrainLabels;
        train_unique = tl_s.umap_train_idx;   % (1 x n_unique) logical
        test_unique  = ~train_unique;          % (1 x n_unique) logical
        test_mask = test_unique(nf.all_to_unique);  % (1 x N_full) logical
    end

% ── Helper: test-unit mask in unique-unit space (for UMAP scatter) ────────────
    function test_unique = getTestMaskUnique()
        tl_s        = ctc.TrainLabels;
        test_unique = ~tl_s.umap_train_idx;   % (1 x n_unique) logical
    end

% ── Tile (1,1): Confidence distribution ──────────────────────────────────────
nexttile(1);
try
    conf      = ctc.UnitConfidence;
    test_mask = getTestMask();
    conf_test = conf(test_mask);
    conf_test = conf_test(isfinite(conf_test));

    if isempty(conf_test)
        error('No test-unit confidence values');
    end

    histogram(conf_test, 20, 'BinLimits', [0 1], 'FaceColor', C_EXC, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.7);

    pct_high = 100 * mean(conf_test > 0.8);
    pct_low  = 100 * mean(conf_test < 0.5);
    text(0.05, 0.93, sprintf('High confidence (>0.8): %.0f%%', pct_high), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    text(0.05, 0.82, sprintf('Low confidence (<0.5): %.0f%%', pct_low), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);

    xlabel('kNN confidence');
    ylabel('Count');
    title('Classification confidence (test units)');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (1,2): Per-chip inhibitory fraction ─────────────────────────────────
nexttile(2);
try
    labels   = ctc.UnitLabels;
    unit_tbl = ctc.FeatureStore.UnitTable;

    if ~ismember('ChipID', unit_tbl.Properties.VariableNames)
        error('ChipID column not found in UnitTable');
    end

    chip_ids    = unit_tbl.ChipID;
    unique_chips = unique(chip_ids, 'stable');
    n_chips      = numel(unique_chips);
    fracs        = zeros(1, n_chips);

    for ci = 1:n_chips
        chip_mask    = chip_ids == unique_chips(ci);
        n_total_chip = sum(chip_mask);
        n_inh_chip   = sum(labels(chip_mask) == 2);
        fracs(ci)    = n_inh_chip / max(n_total_chip, 1);
    end

    bar(1:n_chips, fracs, 'FaceColor', C_INH, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    max_y = max(max(fracs) * 1.2, 0.5);
    patch([0.5, n_chips+0.5, n_chips+0.5, 0.5], [0.15 0.15 0.25 0.25], ...
        [0.2 0.8 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    yline(0.2, '--r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    hold off;

    xticks(1:n_chips);
    xticklabels(string(unique_chips));
    xtickangle(45);
    ylim([0, max_y]);
    ylabel('Inhibitory fraction');
    title('I/E fraction per chip');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (1,3): Supervised UMAP test embedding colored by confidence ──────────
nexttile(3);
try
    red_train = ctc.Reduction.Train;
    red_test  = ctc.Reduction.Test;

    if isempty(red_test)
        error('Reduction.Test not available');
    end

    % Use unique-unit confidence (matches Reduction.Test dimensions).
    % UnitConfidence is broadcast from unique→full; de-duplicate via unique_to_rep.
    nf_loc       = ctc.NormalizedFeatures;
    unique_conf  = ctc.UnitConfidence(nf_loc.unique_to_rep);   % (1 x n_unique)
    test_unique  = getTestMaskUnique();                        % (1 x n_unique)
    conf_test    = unique_conf(test_unique);                   % (1 x n_test_unique)

    scatter(red_test(:, 1), red_test(:, 2), 20, conf_test, 'filled', ...
        'MarkerFaceAlpha', 0.7);
    colormap(gca, 'parula');
    cb = colorbar;
    cb.Label.String = 'Confidence';
    clim([0 1]);

    if ~isempty(red_train)
        hold on;
        scatter(red_train(:, 1), red_train(:, 2), 8, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.5, 'DisplayName', 'Training');
        hold off;
    end

    xlabel('UMAP 1');
    ylabel('UMAP 2');
    title('Test confidence (UMAP)');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (2,1): Mean waveform by predicted class ──────────────────────────────
nexttile(4);
try
    wf = ctc.HarmonizedWaveforms;
    sr = ctc.HarmonizedSR;
    if isempty(wf)
        error('HarmonizedWaveforms not populated');
    end

    labels  = ctc.UnitLabels;
    n_samp  = size(wf, 1);
    t_ms    = (0:(n_samp-1)) / sr * 1000;

    wf_exc = wf(:, labels == 1);
    wf_inh = wf(:, labels == 2);

    hold on;
    if ~isempty(wf_exc)
        mu_exc  = mean(wf_exc, 2);
        sem_exc = std(wf_exc, 0, 2) / sqrt(size(wf_exc, 2));
        shadedErrorBar_local(gca, t_ms, mu_exc', sem_exc', C_EXC, 0.25);
        w_exc = peakToTroughWidth(mu_exc, sr);
        plot(t_ms, mu_exc, '-', 'Color', C_EXC, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('Excitatory (%.2f ms)', w_exc));
    end
    if ~isempty(wf_inh)
        mu_inh  = mean(wf_inh, 2);
        sem_inh = std(wf_inh, 0, 2) / sqrt(size(wf_inh, 2));
        shadedErrorBar_local(gca, t_ms, mu_inh', sem_inh', C_INH, 0.25);
        w_inh = peakToTroughWidth(mu_inh, sr);
        plot(t_ms, mu_inh, '-', 'Color', C_INH, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('Inhibitory (%.2f ms)', w_inh));
    end
    hold off;

    xlabel('Time (ms)');
    ylabel('Amplitude (a.u.)');
    title('Mean waveform by class');
    legend('Location', 'best', 'Box', 'off');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Waveforms unavailable:\n%s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
end

% ── Tile (2,2): Mean ACG by predicted class ───────────────────────────────────
nexttile(5);
try
    acg = ctc.HarmonizedACGs;
    if isempty(acg)
        error('HarmonizedACGs not populated');
    end

    labels      = ctc.UnitLabels;
    acg_bin_ms  = ctc.Parameters.Harmonization.ACGBinSize * 1000;  % s -> ms
    n_bins      = size(acg, 1);
    lag_ms      = ((0:(n_bins-1)) - floor(n_bins/2)) * acg_bin_ms;

    acg_exc = acg(:, labels == 1);
    acg_inh = acg(:, labels == 2);

    hold on;
    if ~isempty(acg_exc)
        mu_exc  = mean(acg_exc, 2);
        sem_exc = std(acg_exc, 0, 2) / sqrt(size(acg_exc, 2));
        shadedErrorBar_local(gca, lag_ms, mu_exc', sem_exc', C_EXC, 0.25);
        plot(lag_ms, mu_exc, '-', 'Color', C_EXC, 'LineWidth', 1.8, ...
            'DisplayName', 'Excitatory');
    end
    if ~isempty(acg_inh)
        mu_inh  = mean(acg_inh, 2);
        sem_inh = std(acg_inh, 0, 2) / sqrt(size(acg_inh, 2));
        shadedErrorBar_local(gca, lag_ms, mu_inh', sem_inh', C_INH, 0.25);
        plot(lag_ms, mu_inh, '-', 'Color', C_INH, 'LineWidth', 1.8, ...
            'DisplayName', 'Inhibitory');
    end
    hold off;

    xlabel('Lag (ms)');
    ylabel('Counts (norm.)');
    title('Mean ACG by class');
    legend('Location', 'best', 'Box', 'off');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('ACGs unavailable:\n%s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
end

% ── Tile (2,3): Feature effect sizes — training vs prediction ─────────────────
nexttile(6);
try
    X_feat = buildXFeat();
    labels = ctc.UnitLabels;
    np     = ctc.NormalizationParams;

    % Prediction Cohen's d (unique units)
    nf         = ctc.NormalizedFeatures;
    % Map full labels -> unique labels using representative index
    labels_unique = labels(nf.unique_to_rep);

    mask1_pred = labels_unique == 1;
    mask2_pred = labels_unique == 2;

    [d_pred, feat_order] = cohensD(X_feat, mask1_pred, mask2_pred, 6);

    feat_names = np.feat_names_trimmed;
    fn_top     = feat_names(feat_order);
    for fi = 1:numel(fn_top)
        s = char(fn_top(fi));
        if numel(s) > 20
            fn_top(fi) = string(s(end-19:end));
        end
    end

    % Training Cohen's d (if available)
    d_train = [];
    tl_s    = ctc.TrainLabels;
    if ~isempty(tl_s) && isfield(tl_s, 'sorted_train_ids')
        train_idx = tl_s.sorted_train_ids;
        y_train   = tl_s.sorted_y_train;
        X_train   = X_feat(train_idx, :);
        mask1_tr  = y_train == 1;
        mask2_tr  = y_train == 2;
        d_tr_full = cohensD_full(X_train, mask1_tr, mask2_tr);
        d_train   = d_tr_full(feat_order);
    end

    n_feat = numel(d_pred);
    y_pos  = 1:n_feat;
    bar_h  = 0.35;

    hold on;
    for fi = 1:n_feat
        col_pred = C_INH;
        if d_pred(fi) < 0; col_pred = C_EXC; end
        barh(y_pos(fi) + bar_h/2, d_pred(fi), bar_h, ...
            'FaceColor', col_pred, 'EdgeColor', 'none', ...
            'FaceAlpha', 1.0, 'DisplayName', '');
        if ~isempty(d_train)
            col_tr = C_INH + 0.4 * (1 - C_INH);
            if d_train(fi) < 0; col_tr = C_EXC + 0.4 * (1 - C_EXC); end
            barh(y_pos(fi) - bar_h/2, d_train(fi), bar_h, ...
                'FaceColor', col_tr, 'EdgeColor', 'none', ...
                'FaceAlpha', 0.6, 'DisplayName', '');
        end
    end

    % Legend proxies
    patch(nan, nan, C_INH, 'FaceAlpha', 1.0, 'EdgeColor', 'none', ...
        'DisplayName', 'Prediction');
    if ~isempty(d_train)
        patch(nan, nan, C_INH + 0.4*(1 - C_INH), 'FaceAlpha', 0.6, ...
            'EdgeColor', 'none', 'DisplayName', 'Training');
    end
    hold off;

    xline(0, '--k', 'LineWidth', 1);
    yticks(y_pos);
    yticklabels(fn_top);
    xlabel("Cohen's d (inh - exc)");
    title('Feature effect sizes (train vs pred)');
    legend('Location', 'best', 'Box', 'off', 'FontSize', 8);
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Save ──────────────────────────────────────────────────────────────────────
ctc.saveDiagnosticFigure(fig, 'classification');

end

%% ── Local helpers ─────────────────────────────────────────────────────────────

function shadedErrorBar_local(ax, x, y, err, color, alpha_val)
% Draw a filled patch for mean +/- err on axes ax. x, y, err are 1xN row vectors.
    x       = x(:)';
    y       = y(:)';
    err     = err(:)';
    x_patch = [x, fliplr(x)];
    y_patch = [y + err, fliplr(y - err)];
    patch(ax, x_patch, y_patch, color, ...
        'FaceAlpha', alpha_val, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function width_ms = peakToTroughWidth(wf_col, sr)
% Estimate peak-to-trough width in ms at half-maximum of abs(wf).
    try
        av = abs(wf_col(:));
        half_max = max(av) / 2;
        above    = av >= half_max;
        idx      = find(above);
        if numel(idx) >= 2
            width_ms = (idx(end) - idx(1)) / sr * 1000;
        else
            width_ms = NaN;
        end
    catch
        width_ms = NaN;
    end
end

function [d_top, feat_order] = cohensD(X, mask1, mask2, n_top)
% Return top-n Cohen's d values and their column indices.
    d_all    = cohensD_full(X, mask1, mask2);
    [~, ord] = sort(abs(d_all), 'descend');
    feat_order = ord(1:min(n_top, numel(ord)));
    d_top      = d_all(feat_order);
end

function d = cohensD_full(X, mask1, mask2)
% Cohen's d for each column of X between two groups.
    n1 = sum(mask1);
    n2 = sum(mask2);
    if n1 < 2 || n2 < 2
        d = zeros(1, size(X, 2));
        return
    end
    mu1    = mean(X(mask1, :), 1);
    mu2    = mean(X(mask2, :), 1);
    s1     = std(X(mask1, :), 0, 1);
    s2     = std(X(mask2, :), 0, 1);
    pooled = sqrt(((n1-1).*s1.^2 + (n2-1).*s2.^2) ./ (n1 + n2 - 2));
    pooled(pooled == 0) = 1;
    d = (mu2 - mu1) ./ pooled;
end
