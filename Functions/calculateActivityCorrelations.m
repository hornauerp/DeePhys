function pop_corr = calculateActivityCorrelations(pop_act, unit_act_mat, pop_sizes, unit_pop_idx)
% CALCULATEACTIVITYCORRELATIONS  Compute correlation between unit activity and population traces.
%
% For each unit, computes the Pearson correlation coefficient with each of the
% provided population activity vectors (e.g. mean excitatory or inhibitory
% firing rate traces).
%
% When pop_sizes and unit_pop_idx are provided, applies leave-one-out correction:
% if unit u belongs to population k, the population mean is recomputed excluding u.
% This removes the trivial self-correlation bias that occurs when a unit is
% included in the population mean it is being correlated with.
%
% INPUTS:
%   pop_act      — (1 x K) cell array of (1 x T) population activity row vectors
%                   (unsmoothed mean firing rates per population)
%   unit_act_mat — (N_units x T) spike count matrix
%   pop_sizes    — (optional) (1 x K) number of units in each population trace.
%                   Required for leave-one-out correction.
%   unit_pop_idx — (optional) (N_units x 1) or (1 x N_units) which population each
%                   unit belongs to (1..K); 0 = does not belong to any population.
%                   Required for leave-one-out correction.
%
% OUTPUTS:
%   pop_corr — (1 x K) cell array; each cell is a (1 x N_units) vector of
%              correlation coefficients

n_units = size(unit_act_mat, 1);
n_pops  = numel(pop_act);

do_loo = nargin >= 4 && ~isempty(pop_sizes) && ~isempty(unit_pop_idx);

pop_corr = arrayfun(@(k) nan(1, n_units), 1:n_pops, 'UniformOutput', false);

for u = 1:n_units
    unit_act = unit_act_mat(u,:);
    for k = 1:n_pops
        if do_loo && unit_pop_idx(u) == k && pop_sizes(k) > 1
            % Leave-one-out: subtract this unit's contribution from the population mean
            loo_mean = (pop_sizes(k) * pop_act{k} - unit_act) / (pop_sizes(k) - 1);
            pop_corr{k}(u) = 1 - pdist2(unit_act, loo_mean, 'correlation');
        else
            pop_corr{k}(u) = 1 - pdist2(unit_act, pop_act{k}, 'correlation');
        end
    end
end
end
