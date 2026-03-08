function computeCorrelations(eia)
% COMPUTECORRELATIONS  Correlate each unit's activity with population E and I traces.
%
% For each culture, computes the Pearson correlation of every unit's binned
% firing rate with the mean inhibitory (I) and mean excitatory (E) population
% activity vectors.
%
% Requires eia.Activity and eia.BinnedMats (run computeActivity first).
% Sets: eia.Correlations

arguments
    eia EIAnalyzer
end

assert(~isempty(eia.Activity), 'Run computeActivity() before computeCorrelations()');

ctc  = eia.Classifier;
rg   = ctc.RecordingGroup;
[cluster_idx, ~] = rg.combineMetadataIndices(ctc.UnitList, ["ChipID", "RecordingDate"]);

corrs = repmat(struct('pop_corr',{},'cell_id',{}), 1, length(eia.Activity));

for c = 1:length(eia.Activity)
    act      = eia.Activity(c);
    pop_act  = {act.I, act.E};
    pop_corr = calculateActivityCorrelations(pop_act, eia.BinnedMats{c});

    % Build cell-type ID vector for this culture's binned rows (1=exc, 2=inh)
    n_rows   = size(eia.BinnedMats{c}, 1);
    cell_id  = ones(1, n_rows);
    culture_labels = ctc.UnitLabels(cluster_idx == c);
    if ~isempty(culture_labels) && length(culture_labels) == n_rows
        cell_id(culture_labels == 2) = 2;
    end

    corrs(c).pop_corr = pop_corr;
    corrs(c).cell_id  = cell_id;
end

eia.Correlations = corrs;
end
