function extractBurstCutouts(eia)
% EXTRACTBURSTCUTOUTS  Extract, align, and cluster burst cutouts.
%
% Two-stage pipeline per culture:
%
%   Stage 0 — Detection & initial alignment
%     1. Detect burst peaks in total activity (optionally masked to burst-state bins)
%     2. Extract +/-PeakCutout bins around each peak
%     3. xcorr-align all bursts to the global mean shape (shift clamped to MaxAlignShift)
%
%   Stage 1 — Outlier rejection (deterministic, iterative)
%     4. Iteratively remove bursts whose Pearson r with the mean template
%        falls below CorrelationThreshold; recompute template each iteration.
%        Cultures with fewer than MinBursts surviving are skipped with a warning.
%
%   Stage 2 — Hierarchical clustering of accepted bursts
%     5. Compute pairwise distance matrix (ClusterDistance metric)
%     6. Average-linkage hierarchical clustering
%     7. Try k = 1..MaxClusters; pick k with highest mean silhouette score.
%        Fall back to k=1 if max silhouette < MinSilhouette or too few bursts.
%     8. Per-cluster xcorr re-alignment to each cluster's own mean
%
% Output fields per culture:
%   .all        — struct (total, i_mat, e_mat): all bursts after step 3
%   .accepted   — struct: bursts surviving outlier rejection (stage 1)
%   .rejected   — struct: discarded outliers
%   .clusters   — (1 x K) struct array: per-cluster aligned bursts
%   .cluster_idx — (1 x N_accepted) cluster assignment for each accepted burst
%   .silhouette — mean silhouette score for the chosen k
%
% Clears NormalizedCutouts to prevent stale results.
%
% Requires eia.Activity and eia.BurstState (run computeActivity + detectBursts first).
% Sets: eia.BurstCutouts

arguments
    eia EIAnalyzer
end

assert(~isempty(eia.BurstState), 'Run detectBursts() before extractBurstCutouts()');

% Invalidate downstream
eia.NormalizedCutouts = [];

p = eia.Parameters.BurstDetection;

empty_cutout = struct( ...
    'all',       struct('total',[],'i_mat',[],'e_mat',[]), ...
    'accepted',  struct('total',[],'i_mat',[],'e_mat',[]), ...
    'rejected',  struct('total',[],'i_mat',[],'e_mat',[]), ...
    'clusters',  struct('total',[],'i_mat',[],'e_mat',[]), ...
    'cluster_idx', [], ...
    'silhouette',  NaN);

cutouts = repmat(empty_cutout, 1, numel(eia.Activity));

for c = 1:numel(eia.Activity)
    act = eia.Activity(c);
    z   = eia.BurstState{c};

    % ── Stage 0: Peak detection ───────────────────────────────────────────
    if p.MaskPeaks
        peak_trace = act.total .* (z == 2);
    else
        peak_trace = act.total;
    end

    [~, locs] = findpeaks(peak_trace, ...
        'MinPeakHeight',     p.PeakHeight, ...
        'MinPeakProminence', p.PeakProminence, ...
        'MinPeakDistance',   p.MinPeakDistance);

    % Remove peaks too close to the edges
    locs(locs <= p.PeakCutout | locs > numel(peak_trace) - p.PeakCutout) = [];

    if numel(locs) < p.MinBursts
        warning('EIAnalyzer:tooFewBursts', ...
            'Culture %d: only %d burst peaks detected (MinBursts=%d) — skipping.', ...
            c, numel(locs), p.MinBursts);
        continue
    end

    % ── Stage 0: Extract cutouts ──────────────────────────────────────────
    n_bursts  = numel(locs);
    w         = 2 * p.PeakCutout + 1;
    burst_mat = nan(n_bursts, w);
    i_mat     = nan(n_bursts, w);
    e_mat     = nan(n_bursts, w);

    for l = 1:n_bursts
        idx            = locs(l) - p.PeakCutout : locs(l) + p.PeakCutout;
        burst_mat(l,:) = act.total(idx);
        i_mat(l,:)     = act.I(idx);
        e_mat(l,:)     = act.E(idx);
    end

    % Drop bursts with NaN bins
    valid     = ~any(isnan(burst_mat), 2);
    burst_mat = burst_mat(valid,:);
    i_mat     = i_mat(valid,:);
    e_mat     = e_mat(valid,:);

    if size(burst_mat, 1) < p.MinBursts
        warning('EIAnalyzer:tooFewBursts', ...
            'Culture %d: only %d valid bursts after NaN removal (MinBursts=%d) — skipping.', ...
            c, size(burst_mat, 1), p.MinBursts);
        continue
    end

    % ── Stage 0: xcorr alignment to global mean (clamped shift) ──────────
    [burst_mat, i_mat, e_mat] = xcorrAlign(burst_mat, i_mat, e_mat, p.PeakCutout, p.MaxAlignShift);

    cutouts(c).all = struct('total', burst_mat, 'i_mat', i_mat, 'e_mat', e_mat);

    % ── Stage 1: Outlier rejection ────────────────────────────────────────
    keep = rejectOutliers(burst_mat, p);

    if sum(keep) < p.MinBursts
        warning('EIAnalyzer:tooFewBursts', ...
            'Culture %d: only %d bursts survive outlier rejection (MinBursts=%d) — keeping top %d.', ...
            c, sum(keep), p.MinBursts, p.MinBursts);
        % Relax: keep the top-MinBursts most correlated bursts
        template = mean(burst_mat, 1, 'omitnan');
        r = zeros(size(burst_mat, 1), 1);
        for b = 1:size(burst_mat, 1)
            row = burst_mat(b,:);
            valid_bins = ~isnan(row) & ~isnan(template);
            if sum(valid_bins) >= 10
                cc = corrcoef(row(valid_bins), template(valid_bins));
                r(b) = cc(1,2);
            end
        end
        [~, rank] = sort(r, 'descend');
        keep = false(size(keep));
        keep(rank(1:min(p.MinBursts, end))) = true;
    end

    reject = ~keep;
    cutouts(c).accepted = struct( ...
        'total', burst_mat(keep,:), 'i_mat', i_mat(keep,:), 'e_mat', e_mat(keep,:));
    cutouts(c).rejected = struct( ...
        'total', burst_mat(reject,:), 'i_mat', i_mat(reject,:), 'e_mat', e_mat(reject,:));

    % ── Stage 2: Hierarchical clustering ─────────────────────────────────
    acc_mat = burst_mat(keep,:);
    acc_i   = i_mat(keep,:);
    acc_e   = e_mat(keep,:);

    % Compute temporal features for clustering (optional, off by default)
    temporal_feats = [];
    if p.TemporalFeatures
        kept_locs = locs(valid);  % locs after NaN removal
        kept_locs = kept_locs(keep);  % locs after outlier rejection
        temporal_feats = buildTemporalFeatures(kept_locs, acc_mat, eia.Parameters.Activity.BinSize);
    end

    [clusters, cluster_idx, sil] = clusterBursts(acc_mat, acc_i, acc_e, p, temporal_feats);

    cutouts(c).clusters    = clusters;
    cutouts(c).cluster_idx = cluster_idx;
    cutouts(c).silhouette  = sil;
end

eia.BurstCutouts = cutouts;
end

%% ── Helper functions ────────────────────────────────────────────────────────

function [burst_mat, i_mat, e_mat] = xcorrAlign(burst_mat, i_mat, e_mat, peak_cutout, max_shift)
% Align bursts to population mean via cross-correlation with clamped shift.
    norm_burst = burst_mat ./ max(burst_mat, [], 2);
    mean_burst = mean(norm_burst, 1, 'omitnan');
    for b = 1:size(burst_mat, 1)
        xc           = xcorr(norm_burst(b,:), mean_burst, peak_cutout);
        [~, max_lag] = max(xc);
        raw_shift    = (peak_cutout + 1) - max_lag;
        shift        = max(-max_shift, min(max_shift, raw_shift));
        burst_mat(b,:) = nanShift(burst_mat(b,:), shift);
        i_mat(b,:)     = nanShift(i_mat(b,:),     shift);
        e_mat(b,:)     = nanShift(e_mat(b,:),     shift);
    end
end

function keep = rejectOutliers(burst_mat, p)
% Iterative correlation-based outlier rejection.
% Returns logical keep vector (same size as burst_mat rows).
    keep = true(size(burst_mat, 1), 1);

    for iter = 1:p.MaxFilterIter
        template      = mean(burst_mat(keep,:), 1, 'omitnan');
        norm_template = template ./ max(abs(template));
        if all(norm_template == 0); break; end

        r = nan(size(burst_mat, 1), 1);
        for b = 1:size(burst_mat, 1)
            row        = burst_mat(b,:);
            valid_bins = ~isnan(row) & ~isnan(template);
            if sum(valid_bins) < 10; r(b) = 0; continue; end
            cc   = corrcoef(row(valid_bins), template(valid_bins));
            r(b) = cc(1,2);
        end

        new_keep = r >= p.CorrelationThreshold;
        if isequal(new_keep, keep); break; end
        keep = new_keep;
    end
end

function [clusters, cluster_idx, best_sil] = clusterBursts(acc_mat, acc_i, acc_e, p, temporal_feats)
% Hierarchical clustering of accepted bursts, with silhouette-based k selection.
%
% Returns:
%   clusters    — (1 x k) struct array, each xcorr-aligned to its cluster mean
%   cluster_idx — (1 x N_accepted) assignment vector
%   best_sil    — mean silhouette for chosen k
%
% Note: distance computation uses a NaN-trimmed version of acc_mat (columns
% that are NaN in any burst are excluded). The full NaN-padded matrices are
% stored in the cluster output structs.
%
% temporal_feats: (optional) (N x F_t) z-scored temporal features appended to
% the shape matrix for distance computation only. Not stored in cluster output.

    if nargin < 5; temporal_feats = []; end

    n = size(acc_mat, 1);

    % Need at least 2 bursts for any clustering; otherwise trivially k=1
    if n < 2 || p.MaxClusters <= 1
        clusters    = struct('total', acc_mat, 'i_mat', acc_i, 'e_mat', acc_e);
        cluster_idx = ones(1, n);
        best_sil    = NaN;
        return
    end

    % Trim to NaN-free columns for distance computation.
    % nanShift pads shifted rows with NaN at the edges; pdist/silhouette
    % return NaN distances for any row containing NaN, corrupting clustering.
    valid_cols  = ~any(isnan(acc_mat), 1);
    acc_mat_trim = acc_mat(:, valid_cols);

    % Append temporal features if provided
    if ~isempty(temporal_feats)
        acc_mat_trim = [acc_mat_trim, temporal_feats];
    end

    if size(acc_mat_trim, 2) < 2
        % Too few valid columns to compute meaningful distances
        clusters    = struct('total', acc_mat, 'i_mat', acc_i, 'e_mat', acc_e);
        cluster_idx = ones(1, n);
        best_sil    = NaN;
        return
    end

    % Pairwise distance + hierarchical linkage (deterministic)
    try
        D = pdist(acc_mat_trim, p.ClusterDistance);
        Z = linkage(D, 'average');
    catch ME
        warning('EIAnalyzer:clusterBursts', ...
            'Clustering failed (%s) — using k=1.', ME.message);
        clusters    = struct('total', acc_mat, 'i_mat', acc_i, 'e_mat', acc_e);
        cluster_idx = ones(1, n);
        best_sil    = NaN;
        return
    end

    % Try k = 2..MaxClusters, pick best mean silhouette
    max_k = min(p.MaxClusters, n - 1);  % can't have more clusters than bursts-1

    if max_k < 2
        % Fewer than 3 bursts — no meaningful multi-cluster comparison possible
        clusters    = struct('total', acc_mat, 'i_mat', acc_i, 'e_mat', acc_e);
        cluster_idx = ones(1, n);
        best_sil    = NaN;
        return
    end

    sil_vals = nan(1, max_k);

    for k = 2:max_k
        idx_k = cluster(Z, 'maxclust', k);
        if numel(unique(idx_k)) < 2; continue; end  % degenerate
        try
            s = silhouette(acc_mat_trim, idx_k, p.ClusterDistance);
            sil_vals(k) = mean(s, 'omitnan');
        catch
            % silhouette may fail for edge cases
        end
    end

    [best_sil_k2plus, best_k_offset] = max(sil_vals(2:end), [], 'omitnan');
    best_k = best_k_offset + 1;  % offset because we skipped k=1

    % Accept k>1 only if silhouette is meaningful
    if isnan(best_sil_k2plus) || best_sil_k2plus < p.MinSilhouette
        best_k   = 1;
        best_sil = NaN;
    else
        best_sil = best_sil_k2plus;
    end

    if best_k == 1
        cluster_idx = ones(1, n);
        clusters    = struct('total', acc_mat, 'i_mat', acc_i, 'e_mat', acc_e);
        return
    end

    cluster_idx = cluster(Z, 'maxclust', best_k)';  % (1 x n)

    % Per-cluster xcorr re-alignment to each cluster's own mean
    clusters(best_k) = struct('total', [], 'i_mat', [], 'e_mat', []);
    for k = 1:best_k
        mask   = cluster_idx == k;
        cm     = acc_mat(mask,:);
        ci_mat = acc_i(mask,:);
        ce_mat = acc_e(mask,:);

        if sum(mask) > 1
            % Re-align within cluster (modest shift allowed — half the global max)
            max_shift_local = max(1, floor(p.MaxAlignShift / 2));
            [cm, ci_mat, ce_mat] = xcorrAlignToMean(cm, ci_mat, ce_mat, max_shift_local);
        end

        clusters(k) = struct('total', cm, 'i_mat', ci_mat, 'e_mat', ce_mat);
    end
end

function [burst_mat, i_mat, e_mat] = xcorrAlignToMean(burst_mat, i_mat, e_mat, max_shift)
% xcorr-align rows to the current mean (used for per-cluster realignment).
    mean_burst = mean(burst_mat ./ max(burst_mat, [], 2), 1, 'omitnan');
    n          = size(burst_mat, 1);
    for b = 1:n
        norm_row     = burst_mat(b,:) ./ max(burst_mat(b,:));
        xc           = xcorr(norm_row, mean_burst, max_shift);
        [~, max_lag] = max(xc);
        shift        = (max_shift + 1) - max_lag;
        shift        = max(-max_shift, min(max_shift, shift));
        burst_mat(b,:) = nanShift(burst_mat(b,:), shift);
        i_mat(b,:)     = nanShift(i_mat(b,:),     shift);
        e_mat(b,:)     = nanShift(e_mat(b,:),     shift);
    end
end

function out = nanShift(row, shift)
% NANSHIFT  Shift a row vector by shift bins, NaN-padding the vacated edge.
    n   = numel(row);
    out = nan(1, n);
    if shift == 0
        out = row;
    elseif shift > 0
        out(shift+1:end) = row(1:end-shift);
    else
        out(1:end+shift) = row(1-shift:end);
    end
end

function feats = buildTemporalFeatures(locs, burst_mat, bin_size)
% BUILDTEMPORALFEATURES  Compute temporal context features for burst clustering.
%
% Returns z-scored features (N_bursts x 3):
%   1. Time since previous burst (seconds; NaN for first burst → replaced by median)
%   2. Absolute time in recording (seconds)
%   3. Burst width at half-maximum (bins)
%
% All features are z-scored so they have comparable scale to the burst shape
% features when appended for pdist computation.

    n = numel(locs);

    % Inter-burst interval
    dt_prev = [NaN; diff(locs(:))] * bin_size;  % seconds
    dt_prev(1) = median(dt_prev(2:end), 'omitnan');  % replace first NaN with median IBI

    % Absolute time in recording
    t_abs = locs(:) * bin_size;

    % Width at half-maximum of total activity cutout
    burst_width = zeros(n, 1);
    for b = 1:n
        row = burst_mat(b,:);
        half_max = max(row, [], 'omitnan') / 2;
        above = row >= half_max;
        if any(above)
            burst_width(b) = sum(above);
        end
    end

    feats = [dt_prev, t_abs, burst_width];

    % Z-score each feature (avoid NaN propagation)
    for f = 1:size(feats, 2)
        col = feats(:, f);
        mu  = mean(col, 'omitnan');
        sd  = std(col, 0, 'omitnan');
        if sd > 0
            feats(:, f) = (col - mu) / sd;
        else
            feats(:, f) = 0;
        end
    end
end
