function computeCorrelations(eia)
% COMPUTECORRELATIONS  Correlate each unit's activity with population E and I traces.
%
% For each culture, computes the Pearson correlation of every unit's binned
% firing rate with the mean inhibitory (I) and mean excitatory (E) population
% activity vectors.
%
% Requires eia.Activity (run computeActivity first).
% Sets: eia.Correlations

arguments
    eia EIAnalyzer
end

assert(~isempty(eia.Activity), 'Run computeActivity() before computeCorrelations()');

ctc   = eia.Classifier;
rg    = ctc.RecordingGroup;
p_act = eia.Parameters.Activity;

% Use the same culture-to-unit mapping as computeActivity (match by object identity)
unit_list  = ctc.UnitList;
all_labels = ctc.UnitLabels;

corrs = repmat(struct('pop_corr',{},'cell_id',{}), 1, length(eia.Activity));

for c = 1:length(eia.Activity)
    act          = eia.Activity(c);
    pop_act      = {act.I, act.E};
    culture      = rg.Cultures(c);
    binned_mat_c = culture.getBinnedSpikeMat(p_act.BinSize, p_act.SecCutout);
    pop_corr     = calculateActivityCorrelations(pop_act, binned_mat_c);

    % Build cell-type ID vector by matching culture units against UnitList
    baseline_units = culture.Units;
    n_rows   = size(binned_mat_c, 1);
    cell_id  = ones(1, n_rows);
    for u = 1:min(length(baseline_units), n_rows)
        match = arrayfun(@(x) x == baseline_units(u), unit_list);
        if any(match)
            lbl = all_labels(find(match, 1));
            if lbl == 2
                cell_id(u) = 2;
            end
        end
    end

    corrs(c).pop_corr = pop_corr;
    corrs(c).cell_id  = cell_id;
end

eia.Correlations = corrs;
end
