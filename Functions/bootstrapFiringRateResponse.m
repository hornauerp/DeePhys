function response = bootstrapFiringRateResponse(pre_mat, post_mat, n_iter, alpha)
% BOOTSTRAPFIRINGRATERESPONSE  Test for significant firing rate changes using bootstrap permutation.
%
% Compares mean firing rates in pre_mat vs post_mat for each unit using a bootstrap
% permutation approach. For each unit, the empirical mean difference (post - pre) is
% compared to the null distribution obtained by permuting bins between windows.
%
% INPUTS:
%   pre_mat  - (N_units x N_bins) spike count matrix — baseline window
%   post_mat - (N_units x N_bins) spike count matrix — post-treatment window;
%              must have the same number of bins as pre_mat
%   n_iter   - Number of bootstrap permutations (default 1000)
%   alpha    - Two-tailed significance level (default 0.01)
%
% OUTPUTS:
%   response - struct with fields:
%              .increase  - indices of units with significant firing rate increase
%              .decrease  - indices of units with significant firing rate decrease
%              .unchanged - indices of units with no significant change

arguments
    pre_mat  (:,:) double
    post_mat (:,:) double
    n_iter   (1,1) double = 1000
    alpha    (1,1) double = 0.01
end

assert(size(pre_mat, 2) == size(post_mat, 2), ...
    'pre_mat and post_mat must have the same number of time bins');

n_units  = size(pre_mat, 1);
n_bins   = size(pre_mat, 2);
full_mat = [pre_mat, post_mat];

emp_diff      = mean(post_mat, 2) - mean(pre_mat, 2);
bootstrp_diff = nan(n_units, n_iter);

for i = 1:n_iter
    rand_idx = randperm(2 * n_bins);
    bootstrp_diff(:,i) = mean(full_mat(:, rand_idx(1:n_bins)),        2) ...
                       - mean(full_mat(:, rand_idx(n_bins+1:end)), 2);
end

increase  = [];
decrease  = [];
unchanged = [];

for u = 1:n_units
    pd = fitdist(bootstrp_diff(u,:)', 'Normal');
    if emp_diff(u) > icdf(pd, 1 - alpha/2)
        increase  = [increase,  u]; %#ok<AGROW>
    elseif emp_diff(u) < icdf(pd, alpha/2)
        decrease  = [decrease,  u]; %#ok<AGROW>
    else
        unchanged = [unchanged, u]; %#ok<AGROW>
    end
end

response = struct('increase', increase, 'decrease', decrease, 'unchanged', unchanged);
end
