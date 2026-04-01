function computeActivity(eia)
% COMPUTEACTIVITY  Build E, I, and total firing rate traces for each culture.
%
% Groups units by culture (via FeatureStore.MetadataTable), bins spike trains
% from UnitData.SpikeTimes, separates by cell type label, and computes mean
% firing rates for each population.
%
% Classified units only (label 1 or 2) contribute to total_fr, so that
% total_fr = mean(E_units + I_units) and the I/E ratio is internally consistent.
% The raw per-unit binned matrix is stored in Activity.binned_mat for reuse by
% computeCorrelations (avoids redundant spike binning).
%
% Clears BurstState, BurstCutouts, NormalizedCutouts, and Correlations to
% prevent stale downstream results after re-running with new parameters.
%
% Requires eia.Classifier.UnitLabels to be set.
% Sets: eia.Activity

arguments
    eia EIAnalyzer
end

% Invalidate all downstream results
eia.BurstState        = [];
eia.BurstCutouts      = [];
eia.NormalizedCutouts = [];
eia.Correlations      = [];

ctc = eia.Classifier;
p   = eia.Parameters.Activity;
map = eia.getCultureUnitMapping();

n_cultures = numel(map.unique_cultures);
ud         = ctc.UnitDataArray;

activity = repmat(struct( ...
    'I', [], 'E', [], 'total', [], 'ratio', [], 'x', [], ...
    'I_raw', [], 'E_raw', [], 'binned_mat', [], 'unit_rows', []), ...
    1, n_cultures);

for c = 1:n_cultures
    cid        = map.unique_cultures(c);
    unit_mask  = map.culture_ids == cid;
    unit_rows  = find(unit_mask);
    ud_indices = map.ud_order(unit_rows);
    valid      = ud_indices > 0;
    unit_rows  = unit_rows(valid);
    ud_indices = ud_indices(valid);

    if isempty(ud_indices)
        continue
    end

    % Build binned spike matrix (N_units x N_bins)
    binned_mat = CellTypeClassifier.binnedSpikeMatrix(ud(ud_indices), p.SecCutout, p.BinSize);
    n_bins     = size(binned_mat, 2);
    x          = p.SecCutout(1) + (0:n_bins-1) * p.BinSize;

    mapped_labels = map.labels(unit_rows);
    inh_rows = mapped_labels == 2;
    exc_rows = mapped_labels == 1;
    cls_rows = inh_rows | exc_rows;  % classified units only

    switch p.RateType
        case "mean"
            rate_fn = @(mat) mean(mat, 1);
        case "sum"
            rate_fn = @(mat) sum(mat, 1);
        otherwise
            error('EIAnalyzer:invalidRateType', 'RateType must be "mean" or "sum", got "%s".', p.RateType);
    end

    if any(inh_rows); I_raw = rate_fn(binned_mat(inh_rows, :)); else; I_raw = zeros(1, n_bins); end
    if any(exc_rows); E_raw = rate_fn(binned_mat(exc_rows, :)); else; E_raw = zeros(1, n_bins); end
    % total uses classified units only so that total is consistent with I and E
    if any(cls_rows); total_raw = rate_fn(binned_mat(cls_rows, :)); else; total_raw = zeros(1, n_bins); end

    % I/E ratio — cap Inf at the finite maximum, replace NaN (0/0) with 0
    ratio        = I_raw ./ E_raw;
    finite_ratio = ratio(isfinite(ratio));
    max_ratio    = max(finite_ratio, [], 'omitnan');
    if isempty(max_ratio); max_ratio = 0; end
    ratio(isinf(ratio)) = max_ratio;
    ratio(isnan(ratio)) = 0;

    % Smooth
    sw       = p.SmoothWindow;
    I_fr     = movmean(I_raw,    sw);
    E_fr     = movmean(E_raw,    sw);
    total_fr = movmean(total_raw, sw);
    ratio    = movmean(ratio,    sw);

    activity(c) = struct( ...
        'I',          I_fr, ...
        'E',          E_fr, ...
        'total',      total_fr, ...
        'ratio',      ratio, ...
        'x',          x, ...
        'I_raw',      I_raw, ...
        'E_raw',      E_raw, ...
        'binned_mat', binned_mat, ...
        'unit_rows',  unit_rows);
end

eia.Activity = activity;
end
