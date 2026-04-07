function diagnosticResponsiveUnits(ctc)
% DIAGNOSTICRESPONSIVEUNITS  2x3 diagnostic figure for identifyResponsiveUnits.
%
% Tiles:
%   (a) [1,1] Per-culture responsive fraction bar chart
%   (b) [1,2] Responsiveness strength distribution (responsive vs non-responsive)
%   (c) [1,3] Dose-response curves (FR vs dose per unit)
%   (d) [2,1] Mean waveform +/- SEM by responsiveness
%   (e) [2,2] Mean ACG +/- SEM by responsiveness
%   (f) [2,3] Summary text
%
% Each tile is guarded by try/catch so a single failure does not abort the figure.
% Calls ctc.saveDiagnosticFigure at the end.

arguments
    ctc CellTypeClassifier
end

C_RESP    = [0.8 0.1 0.1];   % responsive: red
C_NONRESP = [0.5 0.5 0.5];   % non-responsive: grey

fig = figure('Visible', 'on');
set(fig, 'Position', [100 100 1200 700]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Responsive units diagnostics', 'FontWeight', 'bold');

resp_idx = ctc.ResponsiveUnitIdx;   % logical 1xN (full array)
strength = ctc.ResponsiveStrength;  % double 1xN

% ── Tile (a): Per-culture responsive fraction ─────────────────────────────────
nexttile(1);
try
    fs           = ctc.FeatureStore;
    meta         = fs.MetadataTable;
    unit_rec_ids = fs.UnitTable.RecordingID;
    culture_ids  = FeatureStore.getCultureIDsForUnits( ...
        unit_rec_ids, meta, ctc.Parameters.CultureKeys);

    unique_cultures = unique(culture_ids, 'stable');
    n_cult          = numel(unique_cultures);
    fractions       = zeros(1, n_cult);
    for ci = 1:n_cult
        mask        = culture_ids == unique_cultures(ci);
        fractions(ci) = sum(resp_idx(mask)) / max(sum(mask), 1);
    end

    [sorted_frac, sort_ord] = sort(fractions);
    sorted_labels = unique_cultures(sort_ord);

    bar(1:n_cult, sorted_frac, 'FaceColor', C_RESP, 'EdgeColor', 'none');
    overall_frac = sum(resp_idx) / max(numel(resp_idx), 1);
    yline(overall_frac, '--r', 'LineWidth', 1.5);
    ylim([0, max(max(sorted_frac) * 1.2, 0.5)]);
    xticks(1:n_cult);
    xticklabels(string(sorted_labels));
    xtickangle(45);
    ylabel('Responsive fraction');
    title('Responsive fraction per culture');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (b): Strength distribution ──────────────────────────────────────────
nexttile(2);
try
    s_resp    = strength(resp_idx);
    s_nonresp = strength(~resp_idx);

    all_vals = strength(isfinite(strength));
    if isempty(all_vals)
        error('No finite strength values');
    end
    edges = linspace(min(all_vals), max(all_vals) + eps, 31);

    hold on;
    if ~isempty(s_nonresp)
        histogram(s_nonresp, edges, 'FaceColor', C_NONRESP, 'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'DisplayName', 'Non-responsive');
    end
    if ~isempty(s_resp)
        histogram(s_resp, edges, 'FaceColor', C_RESP, 'FaceAlpha', 0.5, ...
            'EdgeColor', 'none', 'DisplayName', 'Responsive');
        xline(median(s_resp), '--', 'Color', C_RESP, 'LineWidth', 1.5);
    end
    hold off;
    xlabel('Responsiveness strength (|rho| or effect size)');
    ylabel('Count');
    title('Strength distribution');
    legend('Location', 'northeast', 'Box', 'off');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (c): Dose-response curves ───────────────────────────────────────────
nexttile(3);
try
    detail = ctc.ResponsivenessDetail;
    if isempty(detail) || isempty(detail.fr_matrix) || isempty(detail.dose_values)
        error('No dose-response data');
    end
    fr_mat    = detail.fr_matrix;   % (n_unique x n_doses)
    dose_vals = detail.dose_values; % (1 x n_doses)

    if numel(unique(dose_vals)) < 2
        error('Single dose — no curve to plot');
    end

    % Map responsive flag from full array to unique units via unique_to_rep
    nf = ctc.NormalizedFeatures;
    if ~isempty(nf)
        rep_rows  = nf.unique_to_rep;           % (1 x n_unique) indices into full array
        resp_u    = resp_idx(rep_rows);          % logical (1 x n_unique)
    else
        % Fallback: assume fr_matrix rows match full array order 1:n_unique
        n_u    = size(fr_mat, 1);
        resp_u = resp_idx(1:n_u);
    end

    % Determine x-axis scale
    pos_doses  = dose_vals(dose_vals > 0);
    use_log    = ~isempty(pos_doses) && ...
                 (max(pos_doses) > 10 * min(pos_doses));

    x_vals = dose_vals;

    hold on;
    % Non-responsive units — thin grey
    nr_rows = find(~resp_u);
    for u = nr_rows
        y = fr_mat(u, :);
        plot(x_vals, y, '-', 'Color', [C_NONRESP, 0.2], 'LineWidth', 0.5);
    end
    % Responsive units — red with alpha
    r_rows = find(resp_u);
    for u = r_rows
        y = fr_mat(u, :);
        plot(x_vals, y, '-', 'Color', [C_RESP, 0.5], 'LineWidth', 0.8);
    end
    % Mean line for responsive units
    if ~isempty(r_rows)
        mean_resp = nanmean(fr_mat(r_rows, :), 1);
        plot(x_vals, mean_resp, '-', 'Color', [0.5 0 0], 'LineWidth', 2.5);
    end
    hold off;

    if use_log
        set(gca, 'XScale', 'log');
    end
    xlabel(sprintf('Dose (%s)', ctc.Parameters.UMAP.GroupingVar));
    ylabel('Firing rate (Hz)');
    title('Dose-response curves');
    box off;
catch ME
    cla; axis off;
    if contains(ME.message, 'dose') || contains(ME.message, 'No ')
        text(0.5, 0.5, 'No dose-response data available', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 9);
    else
        text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', 8, 'Color', [0.6 0 0]);
    end
end

% ── Tile (d): Mean waveform by responsiveness ─────────────────────────────────
nexttile(4);
try
    wf = ctc.HarmonizedWaveforms;   % (N_samples x N_units), full array
    sr = ctc.HarmonizedSR;

    if isempty(wf)
        error('HarmonizedWaveforms not populated');
    end

    n_samp = size(wf, 1);
    t_ms   = (0:(n_samp-1)) / sr * 1000;

    wf_resp    = wf(:, resp_idx);
    wf_nonresp = wf(:, ~resp_idx);

    hold on;
    if ~isempty(wf_nonresp)
        mu_nr  = mean(wf_nonresp, 2);
        sem_nr = std(wf_nonresp, 0, 2) / sqrt(size(wf_nonresp, 2));
        shadedErrorBar_local(t_ms, mu_nr', sem_nr', C_NONRESP, 0.3);
        plot(t_ms, mu_nr, '-', 'Color', C_NONRESP, 'LineWidth', 1.5, ...
            'DisplayName', 'Non-responsive');
    end
    if ~isempty(wf_resp)
        mu_r  = mean(wf_resp, 2);
        sem_r = std(wf_resp, 0, 2) / sqrt(size(wf_resp, 2));
        shadedErrorBar_local(t_ms, mu_r', sem_r', C_RESP, 0.3);
        plot(t_ms, mu_r, '-', 'Color', C_RESP, 'LineWidth', 1.5, ...
            'DisplayName', 'Responsive');
    end
    hold off;

    xlabel('Time (ms)');
    ylabel('Amplitude (a.u.)');
    title('Mean waveform');
    legend('Location', 'best', 'Box', 'off');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Waveforms unavailable:\n%s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
end

% ── Tile (e): Mean ACG by responsiveness ──────────────────────────────────────
nexttile(5);
try
    acg = ctc.HarmonizedACGs;   % (N_bins x N_units), full array

    if isempty(acg)
        error('HarmonizedACGs not populated');
    end

    acg_bin_ms = ctc.Parameters.Harmonization.ACGBinSize * 1000;  % s -> ms
    n_bins     = size(acg, 1);
    lag_ms     = ((0:(n_bins-1)) - floor(n_bins/2)) * acg_bin_ms;

    acg_resp    = acg(:, resp_idx);
    acg_nonresp = acg(:, ~resp_idx);

    hold on;
    if ~isempty(acg_nonresp)
        mu_nr  = mean(acg_nonresp, 2);
        sem_nr = std(acg_nonresp, 0, 2) / sqrt(size(acg_nonresp, 2));
        shadedErrorBar_local(lag_ms, mu_nr', sem_nr', C_NONRESP, 0.3);
        plot(lag_ms, mu_nr, '-', 'Color', C_NONRESP, 'LineWidth', 1.5, ...
            'DisplayName', 'Non-responsive');
    end
    if ~isempty(acg_resp)
        mu_r  = mean(acg_resp, 2);
        sem_r = std(acg_resp, 0, 2) / sqrt(size(acg_resp, 2));
        shadedErrorBar_local(lag_ms, mu_r', sem_r', C_RESP, 0.3);
        plot(lag_ms, mu_r, '-', 'Color', C_RESP, 'LineWidth', 1.5, ...
            'DisplayName', 'Responsive');
    end
    hold off;

    xlabel('Lag (ms)');
    ylabel('Counts (norm.)');
    title('Mean ACG');
    legend('Location', 'best', 'Box', 'off');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('ACGs unavailable:\n%s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
end

% ── Tile (f): Summary text ────────────────────────────────────────────────────
nexttile(6);
try
    n_total    = numel(resp_idx);
    n_resp     = sum(resp_idx);
    n_nonresp  = n_total - n_resp;
    pct_resp   = 100 * n_resp / max(n_total, 1);
    method     = ctc.Parameters.Bootstrap.GroundTruthMethod;
    alpha_val  = ctc.Parameters.Bootstrap.Alpha;
    alpha_fc   = ctc.Parameters.Bootstrap.FullCurveAlpha;

    if lower(string(method)) == "full_curve"
        thresh_str = sprintf('FullCurveAlpha = %.3g', alpha_fc);
    elseif lower(string(method)) == "metadata"
        thresh_str = sprintf('LabelField = "%s"', ctc.Parameters.Bootstrap.LabelField);
    else
        thresh_str = sprintf('Alpha = %.3g', alpha_val);
    end

    summary_lines = {
        sprintf('Total units:      %d', n_total)
        sprintf('Responsive:       %d', n_resp)
        sprintf('Non-responsive:   %d', n_nonresp)
        sprintf('Responsive %%:     %.1f%%', pct_resp)
        ''
        sprintf('Method:  %s', method)
        thresh_str
    };

    axis off;
    text(0.05, 0.95, strjoin(summary_lines, '\n'), ...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 10, ...
        'FontName', 'Monospaced', ...
        'Interpreter', 'none');
    title('Summary');
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Save / hide ───────────────────────────────────────────────────────────────
ctc.saveDiagnosticFigure(fig, 'responsive_units');

end

%% ── Local helper: shaded error bar without toolbox dependency ────────────────

function shadedErrorBar_local(x, y, err, color, alpha_val)
% Draw a filled patch for mean +/- err. x, y, err are 1xN row vectors.
    x      = x(:)';
    y      = y(:)';
    err    = err(:)';
    x_patch = [x, fliplr(x)];
    y_patch = [y + err, fliplr(y - err)];
    patch(x_patch, y_patch, color, ...
        'FaceAlpha', alpha_val, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
