function stability = assessStability(ctc, opts)
% ASSESSSTABILITY  Measure label stability across different random seeds.
%
% Runs the full generateTrainLabels + classifyUnits pipeline N times with
% different RNG seeds and measures pairwise agreement between the resulting
% label vectors using the adjusted Rand index (ARI). High ARI (>0.9) indicates
% stable classification; low ARI suggests the hyperparameters are in a fragile
% region of the search space.
%
% Intended as a post-optimization validation step. Run after
% optimizeHyperparameters to verify the selected parameters produce consistent
% results.
%
% USAGE:
%   ctc.identifyResponsiveUnits();
%   results = ctc.optimizeHyperparameters();
%   ctc.Parameters = parseStructParameters(ctc.Parameters, results.bestParams);
%   stability = ctc.assessStability(NRuns=10);
%
% NAME-VALUE:
%   NRuns   - Number of pipeline runs with different seeds (default 5)
%   Verbose - Print progress (default true)
%
% OUTPUTS:
%   stability - struct with fields:
%     .pairwiseARI   - (N x N) matrix of adjusted Rand indices between runs
%     .meanARI       - scalar mean of off-diagonal ARI values
%     .stdARI        - scalar std of off-diagonal ARI values
%     .labelMatrix   - (N x N_units) matrix of labels from each run
%     .seeds         - (1 x N) seeds used
%     .fractionInh   - (1 x N) inhibitory fraction per run
%     .isStable      - logical, true if meanARI > 0.9

arguments
    ctc CellTypeClassifier
    opts.NRuns   (1,1) double {mustBePositive, mustBeInteger} = 5
    opts.Verbose (1,1) logical = true
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before assessStability()');

N = opts.NRuns;
n_units = numel(ctc.UnitDataArray);

% Generate distinct seeds
base_seed = ctc.Parameters.RNGSeed;
seeds = base_seed + (0:N-1) * 1000;

label_matrix  = nan(N, n_units);
fraction_inh  = nan(1, N);

if opts.Verbose
    fprintf('Assessing stability across %d random seeds...\n', N);
end

for r = 1:N
    temp_params = ctc.Parameters;
    temp_params.RNGSeed = seeds(r);

    temp_ctc = CellTypeClassifier(ctc.FeatureStore, ctc.UnitDataArray, temp_params);
    temp_ctc.ResponsiveUnitIdx = ctc.ResponsiveUnitIdx;

    try
        temp_ctc.generateTrainLabels();
        temp_ctc.classifyUnits();
        label_matrix(r, :) = temp_ctc.UnitLabels;

        valid = ~isnan(temp_ctc.UnitLabels);
        fraction_inh(r) = sum(temp_ctc.UnitLabels(valid) == 2) / sum(valid);

        if opts.Verbose
            fprintf('  Run %d/%d (seed=%d): %.1f%% inhibitory\n', ...
                r, N, seeds(r), 100*fraction_inh(r));
        end
    catch ME
        warning('CellTypeClassifier:assessStability', ...
            'Run %d (seed=%d) failed: %s', r, seeds(r), ME.message);
    end
end

% ── Compute pairwise adjusted Rand index ─────────────────────────────────────
ari_matrix = nan(N, N);
for i = 1:N
    for j = i:N
        if i == j
            ari_matrix(i, j) = 1;
        else
            % Only compare units that have labels in both runs
            valid = ~isnan(label_matrix(i, :)) & ~isnan(label_matrix(j, :));
            if sum(valid) < 10
                ari_matrix(i, j) = NaN;
                ari_matrix(j, i) = NaN;
            else
                ari_matrix(i, j) = adjustedRandIndex(label_matrix(i, valid), label_matrix(j, valid));
                ari_matrix(j, i) = ari_matrix(i, j);
            end
        end
    end
end

% Off-diagonal values
mask = triu(true(N), 1);
off_diag = ari_matrix(mask);
off_diag = off_diag(~isnan(off_diag));

stability.pairwiseARI  = ari_matrix;
stability.meanARI      = mean(off_diag);
stability.stdARI       = std(off_diag);
stability.labelMatrix  = label_matrix;
stability.seeds        = seeds;
stability.fractionInh  = fraction_inh;
stability.isStable     = stability.meanARI > 0.9;

if opts.Verbose
    fprintf('\nStability results:\n');
    fprintf('  Mean ARI:     %.3f +/- %.3f\n', stability.meanARI, stability.stdARI);
    fprintf('  Inh fraction: %.1f%% +/- %.1f%%\n', ...
        100*mean(fraction_inh, 'omitnan'), 100*std(fraction_inh, 'omitnan'));
    if stability.isStable
        fprintf('  Classification is STABLE (ARI > 0.9)\n');
    else
        fprintf('  Classification is UNSTABLE (ARI <= 0.9) — consider re-optimizing\n');
    end
end
end


function ari = adjustedRandIndex(labels_a, labels_b)
% ADJUSTEDRANDINDEX  Compute the adjusted Rand index between two clusterings.
%   ARI = 1 indicates perfect agreement, ARI ~ 0 indicates random labeling.
%   Negative values indicate less agreement than expected by chance.

    n = length(labels_a);
    assert(n == length(labels_b), 'Label vectors must have the same length');

    % Build contingency table
    classes_a = unique(labels_a);
    classes_b = unique(labels_b);
    na = length(classes_a);
    nb = length(classes_b);
    contingency = zeros(na, nb);
    for i = 1:na
        for j = 1:nb
            contingency(i, j) = sum(labels_a == classes_a(i) & labels_b == classes_b(j));
        end
    end

    % Compute ARI from contingency table
    sum_comb_c = sum(contingency(:) .* (contingency(:) - 1)) / 2;
    row_sums = sum(contingency, 2);
    col_sums = sum(contingency, 1);
    sum_comb_a = sum(row_sums .* (row_sums - 1)) / 2;
    sum_comb_b = sum(col_sums .* (col_sums - 1)) / 2;
    n_comb = n * (n - 1) / 2;

    expected = sum_comb_a * sum_comb_b / n_comb;
    max_idx  = (sum_comb_a + sum_comb_b) / 2;

    if max_idx == expected
        ari = 1;  % perfect agreement (or degenerate case)
    else
        ari = (sum_comb_c - expected) / (max_idx - expected);
    end
end
