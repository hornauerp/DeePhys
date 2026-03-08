function pop_corr = calculateActivityCorrelations(pop_act, unit_act_mat)
% CALCULATEACTIVITYCORRELATIONS  Compute correlation between unit activity and population traces.
%
% For each unit, computes the Pearson correlation coefficient with each of the
% provided population activity vectors (e.g. mean excitatory or inhibitory
% firing rate traces).
%
% INPUTS:
%   pop_act      - (1 x K) cell array of (1 x T) population activity row vectors
%   unit_act_mat - (N_units x T) spike count matrix
%
% OUTPUTS:
%   pop_corr - (1 x K) cell array; each cell is a (1 x N_units) vector of
%              correlation coefficients

n_units  = size(unit_act_mat, 1);
pop_corr = arrayfun(@(k) nan(1, n_units), 1:length(pop_act), 'UniformOutput', false);

for u = 1:n_units
    unit_act = unit_act_mat(u,:);
    for k = 1:length(pop_corr)
        pop_corr{k}(u) = 1 - pdist2(unit_act, pop_act{k}, 'correlation');
    end
end
end
