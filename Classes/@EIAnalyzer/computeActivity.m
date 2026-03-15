function computeActivity(eia)
% COMPUTEACTIVITY  Build E, I, and total firing rate traces for each culture.
%
% For each culture, bins the spike trains (using Culture.getBinnedSpikeMat),
% separates units by cell type label (1=excitatory, 2=inhibitory), and
% computes mean firing rates for each population, the total, the I/E ratio,
% and a time axis.
%
% Requires eia.Classifier.UnitLabels to be set.
% Sets: eia.Activity

arguments
    eia EIAnalyzer
end

ctc = eia.Classifier;
rg  = ctc.RecordingGroup;
p   = eia.Parameters.Activity;

n_cultures = length(rg.Cultures);

% Build a reliable culture-to-unit mapping by iterating cultures and
% matching their units against ctc.UnitList, rather than relying on
% combineMetadataIndices which may reorder groups.
unit_list = ctc.UnitList;
all_labels = ctc.UnitLabels;

activity = repmat(struct('I',[],'E',[],'total',[],'ratio',[],'x',[]), 1, n_cultures);

for c = 1:n_cultures
    culture    = rg.Cultures(c);
    binned_mat = culture.getBinnedSpikeMat(p.BinSize, p.SecCutout);

    % Find which entries in ctc.UnitList belong to this culture by matching
    % the Unit objects from culture.Recordings(1).Units (baseline units)
    baseline_units = culture.Units;
    baseline_ids   = [baseline_units.unitID];

    % Match each baseline unit against the full UnitList
    mapped_labels = nan(1, length(baseline_ids));
    for u = 1:length(baseline_units)
        match = arrayfun(@(x) x == baseline_units(u), unit_list);
        if any(match)
            mapped_labels(u) = all_labels(find(match, 1));
        end
    end

    inh_rows = mapped_labels == 2;
    exc_rows = mapped_labels == 1;

    n_bins   = size(binned_mat, 2);
    if any(inh_rows); I_fr = mean(binned_mat(inh_rows, :), 1); else; I_fr = zeros(1, n_bins); end
    if any(exc_rows); E_fr = mean(binned_mat(exc_rows, :), 1); else; E_fr = zeros(1, n_bins); end
    total_fr = mean(binned_mat, 1);

    ratio         = I_fr ./ E_fr;
    finite_ratio  = ratio(isfinite(ratio) & ~isnan(ratio));
    max_ratio     = max(finite_ratio);
    if isempty(max_ratio); max_ratio = 0; end   % all-zero E trace → no finite ratio
    ratio(isinf(ratio))  = max_ratio;
    ratio(isnan(ratio))  = 0;

    activity(c).I     = I_fr;
    activity(c).E     = E_fr;
    activity(c).total = total_fr;
    activity(c).ratio = ratio;
    activity(c).x     = (1:size(binned_mat,2)) * p.BinSize;
end

eia.Activity = activity;
fprintf('Computed E/I activity for %i cultures\n', n_cultures);
end
