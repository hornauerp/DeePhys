function computeActivity(eia)
% COMPUTEACTIVITY  Build E, I, and total firing rate traces for each culture.
%
% For each culture, bins the spike trains (using Culture.getBinnedSpikeMat),
% separates units by cell type label (1=excitatory, 2=inhibitory), and
% computes mean firing rates for each population, the total, the I/E ratio,
% and a time axis.
%
% Requires eia.Classifier.UnitLabels to be set.
% Sets: eia.Activity, eia.BinnedMats

arguments
    eia EIAnalyzer
end

ctc = eia.Classifier;
rg  = ctc.RecordingGroup;
p   = eia.Parameters.Activity;

[cluster_idx, ~] = rg.combineMetadataIndices(ctc.UnitList, ["ChipID", "RecordingDate"]);
n_cultures = length(rg.Cultures);

activity    = repmat(struct('I',[],'E',[],'total',[],'ratio',[],'x',[]), 1, n_cultures);
binned_mats = cell(1, n_cultures);

for c = 1:n_cultures
    culture    = rg.Cultures(c);
    binned_mat = culture.getBinnedSpikeMat(p.BinSize, p.SecCutout);
    binned_mats{c} = binned_mat;

    % Map ctc.UnitLabels (indexed into UnitList) onto culture.Units order
    culture_labels = ctc.UnitLabels(cluster_idx == c);
    culture_units  = ctc.UnitList(cluster_idx == c);
    baseline_ids   = [culture.Units.unitID];
    all_ids        = [culture_units.unitID];

    [~, map_idx]   = ismember(baseline_ids, all_ids);
    mapped_labels  = nan(1, length(baseline_ids));
    valid = map_idx > 0;
    mapped_labels(valid) = culture_labels(map_idx(valid));

    inh_rows = mapped_labels == 2;
    exc_rows = mapped_labels == 1;

    I_fr     = mean(binned_mat(inh_rows, :), 1);
    E_fr     = mean(binned_mat(exc_rows, :), 1);
    total_fr = mean(binned_mat, 1);

    if isempty(I_fr) || ~any(inh_rows); I_fr = zeros(1, size(binned_mat, 2)); end
    if isempty(E_fr) || ~any(exc_rows); E_fr = zeros(1, size(binned_mat, 2)); end

    ratio         = I_fr ./ E_fr;
    finite_ratio  = ratio(isfinite(ratio) & ~isnan(ratio));
    max_ratio     = max(finite_ratio);
    ratio(isinf(ratio))  = max_ratio;
    ratio(isnan(ratio))  = 0;

    activity(c).I     = I_fr;
    activity(c).E     = E_fr;
    activity(c).total = total_fr;
    activity(c).ratio = ratio;
    activity(c).x     = (1:size(binned_mat,2)) * p.BinSize;
end

eia.Activity   = activity;
eia.BinnedMats = binned_mats;
fprintf('Computed E/I activity for %i cultures\n', n_cultures);
end
