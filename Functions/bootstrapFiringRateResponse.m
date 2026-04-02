function response = bootstrapFiringRateResponse(pre_mat, post_mat, n_iter, alpha)
% BOOTSTRAPFIRINGRATERESPONSE  Test for significant firing rate changes using bootstrap permutation.
%
% Compares mean firing rates in pre_mat vs post_mat for each unit using a bootstrap
% permutation approach. For each unit, the empirical mean difference (post - pre) is
% compared to a Normal distribution fitted to the bootstrap null distribution.
% This parametric test allows meaningful inference at extreme significance levels
% (e.g. alpha = 1e-10) even with moderate iteration counts (n_iter = 1000).
%
% When n_iter is a struct with MaxNIter and ConvergenceTol fields, the bootstrap
% runs in batches of 500 with early stopping when the Normal fit parameters
% (mu, sigma) stabilize within ConvergenceTol.
%
% INPUTS:
%   pre_mat  - (N_units x N_bins) spike count matrix — baseline window
%   post_mat - (N_units x N_bins) spike count matrix — post-treatment window
%   n_iter   - Number of bootstrap permutations (default 1000), OR a struct with
%              fields .NIter (initial), .MaxNIter (cap), .ConvergenceTol
%   alpha    - Two-tailed significance level (default 1e-10)
%
% OUTPUTS:
%   response - struct with fields:
%              .increase  - indices of units with significant firing rate increase
%              .decrease  - indices of units with significant firing rate decrease
%              .unchanged - indices of units with no significant change
%              .strength  - (1 x N_units) effect size (|emp_diff| / null_std)
%              .p_values  - (1 x N_units) two-tailed p-values from Normal null fit

arguments
    pre_mat  (:,:) double
    post_mat (:,:) double
    n_iter   = 1000
    alpha    (1,1) double = 1e-10
end

assert(size(pre_mat, 2) == size(post_mat, 2), ...
    'pre_mat and post_mat must have the same number of time bins');

n_units  = size(pre_mat, 1);
n_bins   = size(pre_mat, 2);
full_mat = [pre_mat, post_mat];

emp_diff = mean(post_mat, 2) - mean(pre_mat, 2);

% -- Determine iteration scheme ----------------------------------------------
if isstruct(n_iter)
    % Adaptive iterations with convergence monitoring
    init_iter   = n_iter.NIter;
    max_iter    = n_iter.MaxNIter;
    conv_tol    = n_iter.ConvergenceTol;
else
    init_iter   = n_iter;
    max_iter    = n_iter;
    conv_tol    = 0;  % no convergence check when fixed iterations
end

batch_size = 500;

% -- Adaptive bootstrap with convergence monitoring (Phase 3) ---------------
bootstrp_diff = nan(n_units, 0);
total_iter    = 0;
converged     = false;

% First pass: run at least init_iter iterations
n_first = max(init_iter, batch_size);
new_batch = nan(n_units, n_first);
for i = 1:n_first
    rand_idx = randperm(2 * n_bins);
    new_batch(:, i) = mean(full_mat(:, rand_idx(1:n_bins)),    2) ...
                    - mean(full_mat(:, rand_idx(n_bins+1:end)), 2);
end
bootstrp_diff = [bootstrp_diff, new_batch];
total_iter    = total_iter + n_first;

% Additional batches until convergence or max_iter
while total_iter < max_iter && ~converged
    this_batch = min(batch_size, max_iter - total_iter);
    new_batch  = nan(n_units, this_batch);
    for i = 1:this_batch
        rand_idx = randperm(2 * n_bins);
        new_batch(:, i) = mean(full_mat(:, rand_idx(1:n_bins)),    2) ...
                        - mean(full_mat(:, rand_idx(n_bins+1:end)), 2);
    end
    bootstrp_diff = [bootstrp_diff, new_batch]; %#ok<AGROW>
    total_iter    = total_iter + this_batch;

    if conv_tol > 0 && total_iter >= 1000
        % Convergence check: compare full history vs first half
        half     = floor(size(bootstrp_diff, 2) / 2);
        mu_all   = mean(bootstrp_diff, 2);
        sd_all   = std(bootstrp_diff, 0, 2);
        mu_half  = mean(bootstrp_diff(:, 1:half), 2);
        sd_half  = std(bootstrp_diff(:, 1:half), 0, 2);
        ref_mu   = max(abs(mu_all), 1e-10);
        ref_sd   = max(abs(sd_all), 1e-10);
        converged = max(abs(mu_all - mu_half) ./ ref_mu) < conv_tol && ...
                    max(abs(sd_all - sd_half) ./ ref_sd) < conv_tol;
    end
end

if max_iter > init_iter
    if converged
        fprintf('Bootstrap: %d iterations (converged)\n', total_iter);
    else
        fprintf('Bootstrap: %d iterations (max reached)\n', total_iter);
    end
end

% -- Fit Normal null distribution and threshold ------------------------------
increase  = [];
decrease  = [];
unchanged = [];
p_values  = nan(1, n_units);

for u = 1:n_units
    pd    = fitdist(bootstrp_diff(u, :)', 'Normal');
    upper = icdf(pd, 1 - alpha / 2);
    lower = icdf(pd, alpha / 2);

    % Two-tailed p-value for FDR correction support
    p_values(u) = 2 * min(cdf(pd, emp_diff(u)), 1 - cdf(pd, emp_diff(u)));

    if emp_diff(u) > upper
        increase  = [increase,  u]; %#ok<AGROW>
    elseif emp_diff(u) < lower
        decrease  = [decrease,  u]; %#ok<AGROW>
    else
        unchanged = [unchanged, u]; %#ok<AGROW>
    end
end

% Effect size: |empirical difference| normalized by bootstrap null std.
null_std = std(bootstrp_diff, 0, 2);
null_std(null_std == 0) = 1;
strength = (abs(emp_diff) ./ null_std)';

response = struct('increase', increase, 'decrease', decrease, ...
    'unchanged', unchanged, 'strength', strength, 'p_values', p_values);
end
