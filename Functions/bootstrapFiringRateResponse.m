function response = bootstrapFiringRateResponse(pre_mat, post_mat, n_iter, alpha, opts)
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
% NAME-VALUE INPUTS:
%   UnitWisePermutation - logical (default false). When true, each unit gets
%     an independent randperm per iteration. More conservative for correlated
%     spike trains (e.g. synchronous MEA bursting), but ~N_units times slower.
%   NetworkCorrection   - logical (default false). When true, subtract the
%     population-mean FR change before testing. Removes the network-wide
%     disinhibition confound (e.g. GABAergic blockade). Units are then tested
%     for ABOVE-AVERAGE response rather than any increase.
%
% OUTPUTS:
%   response - struct with fields:
%              .increase       - indices of units with significant FR increase
%              .decrease       - indices of units with significant FR decrease
%              .unchanged      - indices of units with no significant change
%              .strength       - (1 x N_units) effect size (|emp_diff| / null_std)
%              .p_values       - (1 x N_units) two-tailed p-values from Normal null
%              .n_non_normal   - number of units with non-Normal null (Jarque-Bera)
%              .frac_non_normal- fraction of units with non-Normal null

arguments
    pre_mat  (:,:) double
    post_mat (:,:) double
    n_iter   = 1000
    alpha    (1,1) double = 1e-10
    opts.UnitWisePermutation (1,1) logical = false
    opts.NetworkCorrection   (1,1) logical = false
end

assert(size(pre_mat, 2) == size(post_mat, 2), ...
    'pre_mat and post_mat must have the same number of time bins');

n_units  = size(pre_mat, 1);
n_bins   = size(pre_mat, 2);
full_mat = [pre_mat, post_mat];

emp_diff = mean(post_mat, 2) - mean(pre_mat, 2);

% -- Network-level firing rate correction ------------------------------------
% When NetworkCorrection=true, subtract the population-mean firing rate change
% before per-unit testing. Units are then tested for above-average response,
% removing the network-wide disinhibition confound (e.g. under GABAergic
% blockade, all units tend to increase FR; only units that increase *more than
% average* are likely to be truly inhibitory).
if opts.NetworkCorrection
    pop_mean_diff = mean(emp_diff);
    emp_diff      = emp_diff - pop_mean_diff;
    % Apply the same correction to the full matrix (post columns only) so
    % that the bootstrap null reflects the corrected distribution.
    full_mat(:, n_bins+1:end) = full_mat(:, n_bins+1:end) - pop_mean_diff;
end

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
if opts.UnitWisePermutation
    % Independent permutation per unit — avoids correlated null distributions
    % for units with shared network-driven variance (e.g. synchronous bursting).
    for i = 1:n_first
        for u = 1:n_units
            rand_idx = randperm(2 * n_bins);
            new_batch(u, i) = mean(full_mat(u, rand_idx(1:n_bins)), 2) ...
                            - mean(full_mat(u, rand_idx(n_bins+1:end)), 2);
        end
    end
else
    % Shared permutation across units (original behavior, faster).
    for i = 1:n_first
        rand_idx = randperm(2 * n_bins);
        new_batch(:, i) = mean(full_mat(:, rand_idx(1:n_bins)),    2) ...
                        - mean(full_mat(:, rand_idx(n_bins+1:end)), 2);
    end
end
bootstrp_diff = [bootstrp_diff, new_batch];
total_iter    = total_iter + n_first;

% Additional batches until convergence or max_iter
while total_iter < max_iter && ~converged
    this_batch = min(batch_size, max_iter - total_iter);
    new_batch  = nan(n_units, this_batch);
    if opts.UnitWisePermutation
        for i = 1:this_batch
            for u = 1:n_units
                rand_idx = randperm(2 * n_bins);
                new_batch(u, i) = mean(full_mat(u, rand_idx(1:n_bins)), 2) ...
                                - mean(full_mat(u, rand_idx(n_bins+1:end)), 2);
            end
        end
    else
        for i = 1:this_batch
            rand_idx = randperm(2 * n_bins);
            new_batch(:, i) = mean(full_mat(:, rand_idx(1:n_bins)),    2) ...
                            - mean(full_mat(:, rand_idx(n_bins+1:end)), 2);
        end
    end
    bootstrp_diff = [bootstrp_diff, new_batch]; %#ok<AGROW>
    total_iter    = total_iter + this_batch;

    if conv_tol > 0 && total_iter >= 1000
        % Convergence check: compare two consecutive non-overlapping halves
        half   = floor(size(bootstrp_diff, 2) / 2);
        mu_h1  = mean(bootstrp_diff(:, 1:half), 2);
        sd_h1  = std(bootstrp_diff(:, 1:half), 0, 2);
        mu_h2  = mean(bootstrp_diff(:, half+1:end), 2);
        sd_h2  = std(bootstrp_diff(:, half+1:end), 0, 2);
        ref_mu = max(max(abs([mu_h1, mu_h2]), [], 2), 1e-10);
        ref_sd = max(max(abs([sd_h1, sd_h2]), [], 2), 1e-10);
        converged = max(abs(mu_h1 - mu_h2) ./ ref_mu) < conv_tol && ...
                    max(abs(sd_h1 - sd_h2) ./ ref_sd) < conv_tol;
    end
end

if max_iter > init_iter
    if converged
        fprintf('Bootstrap: %d iterations (converged)\n', total_iter);
    else
        fprintf('Bootstrap: %d iterations (max reached)\n', total_iter);
    end
end

% -- Gaussianity diagnostic --------------------------------------------------
% The Normal parametric fit extrapolates into the tail at extreme alpha values
% (e.g. 1e-10 ≈ 6.5 sigma). Validate the Gaussian assumption with a
% Jarque-Bera test on each unit's null distribution. When > 20% of units are
% non-Normal AND alpha < 1e-3, fall back to empirical quantiles per unit.
non_normal_mask = false(n_units, 1);
n_non_normal = 0;
for u = 1:n_units
    try
        h = jbtest(bootstrp_diff(u, :)', 0.05);
        if h
            non_normal_mask(u) = true;
            n_non_normal = n_non_normal + 1;
        end
    catch
        % jbtest requires >= 4 samples; skip if not available
    end
end
frac_non_normal = n_non_normal / max(n_units, 1);

use_empirical = non_normal_mask & (frac_non_normal > 0.2) & (alpha < 1e-3);
if any(use_empirical)
    alpha_min_empirical = 1 / total_iter;
    warning('bootstrapFiringRateResponse:nonNormalNull', ...
        ['%.0f%% of units have non-Normal bootstrap null distributions ' ...
         '(Jarque-Bera, alpha=0.05). Using empirical quantiles for those units. ' ...
         'Minimum achievable alpha with %d iterations: %.0e.'], ...
        100 * frac_non_normal, total_iter, alpha_min_empirical);
    if alpha < alpha_min_empirical
        warning('bootstrapFiringRateResponse:empiricalResolution', ...
            ['alpha (%.0e) is below the empirical resolution floor (1/n_iter = %.0e). ' ...
             'Increase n_iter for reliable inference at this significance level.'], ...
            alpha, alpha_min_empirical);
    end
elseif frac_non_normal > 0.2 && alpha < 1e-3
    warning('bootstrapFiringRateResponse:nonNormalNull', ...
        ['%.0f%% of units have non-Normal bootstrap null distributions ' ...
         '(Jarque-Bera, alpha=0.05). Parametric p-values at extreme alpha ' ...
         '(%.0e) may be unreliable. Consider a more moderate alpha or ' ...
         'increasing n_iter.'], 100 * frac_non_normal, alpha);
end

% -- Fit null distribution and threshold per unit ----------------------------
increase  = [];
decrease  = [];
unchanged = [];
p_values  = nan(1, n_units);

for u = 1:n_units
    if use_empirical(u)
        % Empirical quantiles — resolution limited to 1/total_iter
        sorted_null = sort(bootstrp_diff(u, :));
        upper_idx = max(1, round((1 - alpha/2) * total_iter));
        lower_idx = max(1, round(alpha/2 * total_iter));
        upper = sorted_null(min(upper_idx, total_iter));
        lower = sorted_null(lower_idx);
        % Two-tailed empirical p-value
        p_values(u) = 2 * min( ...
            mean(bootstrp_diff(u,:) >= emp_diff(u)), ...
            mean(bootstrp_diff(u,:) <= emp_diff(u)));
    else
        pd    = fitdist(bootstrp_diff(u, :)', 'Normal');
        upper = icdf(pd, 1 - alpha / 2);
        lower = icdf(pd, alpha / 2);
        % Two-tailed p-value for FDR correction support
        p_values(u) = 2 * min(cdf(pd, emp_diff(u)), 1 - cdf(pd, emp_diff(u)));
    end

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
    'unchanged', unchanged, 'strength', strength, 'p_values', p_values, ...
    'n_non_normal', n_non_normal, 'frac_non_normal', frac_non_normal);
end
