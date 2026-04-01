function response = bootstrapDoseResponse(culture, n_iter, alpha)
% BOOTSTRAPDOSERESPONSE  Test for monotonic firing rate changes across the full dose-response curve.
%
% For each unit, computes mean firing rate per recording (concentration level),
% then tests whether firing rate correlates with concentration rank using a
% Spearman rank correlation with bootstrap permutation significance.
%
% This uses the full dose-response trajectory rather than comparing only
% the initial and final segments.
%
% INPUTS:
%   culture - Culture object (must have >=2 recordings with GroupingValues)
%   n_iter  - Number of permutation iterations (default 1000)
%   alpha   - Two-tailed significance level (default 1e-10)
%
% OUTPUTS:
%   response - struct with fields:
%              .increase  - indices of units with significant positive dose-response (inhibitory candidates)
%              .decrease  - indices of units with significant negative dose-response
%              .unchanged - indices of units with no significant trend
%              .strength  - (1 x N_units) absolute Spearman rho for each unit (0-1)
%                           Higher values indicate stronger dose-response monotonicity.

arguments
    culture Culture
    n_iter  (1,1) double = 1000
    alpha   (1,1) double = 1e-10
end

n_recordings = culture.N;
assert(n_recordings >= 2, 'Need at least 2 recordings for dose-response analysis');

baseline_units = culture.Units;
n_units = numel(baseline_units);
baseline_ids = [baseline_units.TemplateID];

% ── Compute mean firing rate per unit per recording ─────────────────────
rates = zeros(n_units, n_recordings);
for r = 1:n_recordings
    rec = culture.Recordings(r);
    rec_ids = [rec.Units.TemplateID];
    [~, row_idx] = ismember(baseline_ids, rec_ids);
    duration = rec.RecordingInfo.Duration;
    for u = 1:n_units
        if row_idx(u) == 0; continue; end
        rates(u, r) = numel(rec.Units(row_idx(u)).SpikeTimes) / duration;
    end
end

% ── Concentration ranks for Spearman correlation ────────────────────────
conc_rank = tiedrank(culture.GroupingValues(:));

% ── Empirical Spearman correlation per unit ─────────────────────────────
emp_corr = zeros(n_units, 1);
for u = 1:n_units
    emp_corr(u) = corr(rates(u,:)', conc_rank, 'type', 'Spearman');
end

% ── Permutation test: shuffle concentration ranks ───────────────────────
boot_corr = zeros(n_units, n_iter);
for i = 1:n_iter
    perm_rank = conc_rank(randperm(n_recordings));
    for u = 1:n_units
        boot_corr(u, i) = corr(rates(u,:)', perm_rank, 'type', 'Spearman');
    end
end

% ── Significance via Normal fit (allows extrapolation at extreme alpha) ─
increase  = [];
decrease  = [];
unchanged = [];

for u = 1:n_units
    pd = fitdist(boot_corr(u,:)', 'Normal');
    upper = icdf(pd, 1 - alpha / 2);
    lower = icdf(pd, alpha / 2);
    if emp_corr(u) > upper
        increase  = [increase,  u]; %#ok<AGROW>
    elseif emp_corr(u) < lower
        decrease  = [decrease,  u]; %#ok<AGROW>
    else
        unchanged = [unchanged, u]; %#ok<AGROW>
    end
end

response = struct( ...
    'increase',  increase, ...
    'decrease',  decrease, ...
    'unchanged', unchanged, ...
    'strength',  abs(emp_corr'));
end
