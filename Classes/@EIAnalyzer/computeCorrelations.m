function computeCorrelations(eia)
% COMPUTECORRELATIONS  Correlate each unit's activity with population E and I traces.
%
% For each culture, computes the Pearson correlation of every unit's binned
% firing rate with the unsmoothed mean inhibitory (I_raw) and excitatory (E_raw)
% population activity vectors stored in Activity.
%
% Reuses Activity.binned_mat so spike trains are not re-binned. Population
% traces used for correlation are the unsmoothed I_raw / E_raw, avoiding the
% asymmetry that would arise from correlating raw unit activity against
% smoothed population means.
%
% Requires eia.Activity (run computeActivity first).
% Sets: eia.Correlations

arguments
    eia EIAnalyzer
end

assert(~isempty(eia.Activity), 'Run computeActivity() before computeCorrelations()');

all_labels = eia.Classifier.UnitLabels;
corrs      = repmat(struct('pop_corr', [], 'cell_id', []), 1, numel(eia.Activity));

for c = 1:numel(eia.Activity)
    fprintf('  computeCorrelations: culture %d/%d\n', c, numel(eia.Activity));
    act = eia.Activity(c);

    if isempty(act.binned_mat)
        continue
    end

    % Build cell-type ID vector (0 = unclassified, 1 = excitatory, 2 = inhibitory)
    % unit_rows indexes into FeatureStore/UnitLabels order — cached from computeActivity
    cell_id = all_labels(act.unit_rows);
    cell_id(isnan(cell_id)) = 0;

    % Compute population MEAN traces directly from binned_mat for correlation.
    % This ensures the leave-one-out formula (which assumes means) is correct
    % regardless of the Activity.RateType setting (mean vs sum).
    inh_mask = cell_id == 2;
    exc_mask = cell_id == 1;
    n_inh    = sum(inh_mask);
    n_exc    = sum(exc_mask);

    if n_inh > 0; I_mean = mean(act.binned_mat(inh_mask, :), 1); else; I_mean = zeros(1, size(act.binned_mat, 2)); end
    if n_exc > 0; E_mean = mean(act.binned_mat(exc_mask, :), 1); else; E_mean = zeros(1, size(act.binned_mat, 2)); end

    % Population membership for leave-one-out correction:
    %   pop_act{1} = I_mean → population index 1 = cell_id 2
    %   pop_act{2} = E_mean → population index 2 = cell_id 1
    pop_idx = zeros(size(cell_id));
    pop_idx(cell_id == 2) = 1;  % inhibitory units belong to I population (pop_act{1})
    pop_idx(cell_id == 1) = 2;  % excitatory units belong to E population (pop_act{2})

    pop_act  = {I_mean, E_mean};
    pop_corr = calculateActivityCorrelations(pop_act, act.binned_mat, [n_inh, n_exc], pop_idx);

    corrs(c).pop_corr = pop_corr;
    corrs(c).cell_id  = cell_id;
end

eia.Correlations = corrs;
end
