function results = optimizeHyperparameters(ctc, opts)
% OPTIMIZEHYPERPARAMETERS  Bayesian optimization of UMAP topology parameters.
%
% Optimizes 3 UMAP topology parameters (MinDist, Spread, SupervisedNDims)
% to produce a well-structured supervised embedding. These directly control
% embedding shape — fixing problems like degenerate ellipses or collapsed
% clusters.
%
% Only topology parameters are optimized. Parameters that affect label
% assignment or class balance (TargetWeight, CounterexampleRatio) are
% excluded because the objective metric can be trivially gamed by these.
%
% Efficiency: training labels and the normalised feature matrix are generated
% once before the optimisation loop. Each iteration only re-runs the
% supervised UMAP (classifyUnits). The high-D kNN for trustworthiness is
% precomputed once since the original feature space is invariant.
%
% Requires identifyResponsiveUnits() to have been run first (ground truth
% is fixed during optimization).
%
% USAGE:
%   ctc = CellTypeClassifier(fs, ud, params);
%   ctc.identifyResponsiveUnits();
%   results = ctc.optimizeHyperparameters();
%   ctc.Parameters = parseStructParameters(ctc.Parameters, results.bestParams);
%   ctc.generateTrainLabels();
%   ctc.classifyUnits();
%
% NAME-VALUE:
%   Verbose - Print progress during optimization (default true)
%
% OBJECTIVE METRIC (set via Parameters.BayesianOptimization.ObjectiveMetric):
%   "trustworthiness" (default) — 1 - T, where T (Venna & Kaski 2006) measures
%       how faithfully the low-D embedding preserves high-D neighborhoods.
%       T = 1 is perfect; works in any output dimension.
%   "silhouette"       — negative mean silhouette on supervised embedding.
%   "combined"         — (1-T) + (-mean_sil).
%
% OPTIMIZED VARIABLES (always 3):
%   MinDist        — UMAP min_dist (log-scaled over MinDistRange)
%   Spread         — UMAP spread (over SpreadRange)
%   SupervisedNDims— supervised UMAP output dimensions (integer, SupervisedNDimsRange)
%
% Fixed: 30 Bayesian evaluations.
%
% OUTPUTS:
%   results - struct with fields:
%     .bestParams     - struct of optimal parameter values (UMAP subfield)
%     .bestObjective  - best (lowest) objective value achieved
%     .bayesoptResult - full BayesianOptimization object
%     .allObjectives  - table of all evaluations
%     .nVars          - number of variables optimized (always 3)

arguments
    ctc CellTypeClassifier
    opts.Verbose (1,1) logical = true
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before optimizeHyperparameters()');

p_bo = ctc.Parameters.BayesianOptimization;

% -- Define 3 optimizable variables ------------------------------------------
vars = [ ...
    optimizableVariable('MinDist',         p_bo.MinDistRange,        'Transform', 'log'), ...
    optimizableVariable('Spread',          p_bo.SpreadRange), ...
    optimizableVariable('SupervisedNDims', p_bo.SupervisedNDimsRange, 'Type', 'integer')];
n_vars    = 3;
max_evals = 30;

% -- Generate training labels once (fixed across all iterations) --------------
% Feature extraction, unsupervised UMAP, outlier detection, and training label
% assembly happen here. Each BayOpt iteration only re-runs the supervised UMAP.
if opts.Verbose
    fprintf('Generating fixed training labels for optimization...\n');
end
ref_ctc = CellTypeClassifier(ctc.FeatureStore, ctc.UnitDataArray, ctc.Parameters);
ref_ctc.ResponsiveUnitIdx       = ctc.ResponsiveUnitIdx;
ref_ctc.ResponsiveUnitDirection = ctc.ResponsiveUnitDirection;
ref_ctc.ResponsiveStrength      = ctc.ResponsiveStrength;
ref_ctc.CounterexampleUnitIdx   = ctc.CounterexampleUnitIdx;
ref_ctc.generateTrainLabels();

fixed_labels = ref_ctc.TrainLabels;
fixed_np     = ref_ctc.NormalizationParams;
fixed_nf     = ref_ctc.NormalizedFeatures;
n_unique     = fixed_labels.n_unique;

% -- Precompute high-D feature matrix and kNN (invariant across iterations) ---
X_all_high = normalize(fixed_nf.X_pergroup, 'center', fixed_np.mu_global, 'scale', fixed_np.sigma_global);
X_all_high(:, fixed_np.nan_cols) = [];
X_all_high = X_all_high ./ fixed_np.scale;
if isfield(fixed_np, 'feature_selection_mask') && ~all(fixed_np.feature_selection_mask)
    X_all_high = X_all_high(:, fixed_np.feature_selection_mask);
end
train_mask = logical(fixed_labels.umap_train_idx);
X_high     = [X_all_high(train_mask, :); X_all_high(~train_mask, :)];

k_trust = 15;
k_hi    = min(5 * k_trust, size(X_high, 1) - 1);
[nn_high_precomputed, ~] = knnsearch(X_high, X_high, 'K', k_hi + 1);
nn_high_precomputed = nn_high_precomputed(:, 2:end);  % exclude self

if opts.Verbose
    fprintf('Precomputed: %d unique units (%d train), %d features, high-D kNN (k_hi=%d)\n', ...
        n_unique, sum(train_mask), size(X_high, 2), k_hi);
end

% -- Objective function -------------------------------------------------------
function loss = objective(x)
    temp_params = ctc.Parameters;
    temp_params.UMAP.MinDist         = x.MinDist;
    temp_params.UMAP.Spread          = x.Spread;
    temp_params.UMAP.SupervisedNDims = x.SupervisedNDims;

    temp_ctc = CellTypeClassifier(ctc.FeatureStore, ctc.UnitDataArray, temp_params);
    temp_ctc.TrainLabels         = fixed_labels;
    temp_ctc.NormalizationParams = fixed_np;
    temp_ctc.NormalizedFeatures  = fixed_nf;
    temp_ctc.Reduction           = struct();  % classifyUnits writes .Train/.Test/.Extras

    try
        temp_ctc.classifyUnits();

        all_reduction = [temp_ctc.Reduction.Train; temp_ctc.Reduction.Test];

        % Reconstruct unique-unit labels (UnitLabels is broadcast to n_units_full;
        % umap_test_idx is n_unique-length, so index via the representative row map)
        unique_labels = temp_ctc.UnitLabels(fixed_nf.unique_to_rep);
        train_labels_vec = unique_labels(train_mask)';
        test_labels_vec  = unique_labels(~train_mask)';
        all_labels       = [train_labels_vec; test_labels_vec];

        valid = ~isnan(all_labels);
        if sum(valid) < 10 || numel(unique(all_labels(valid))) < 2
            loss = 1;
            return
        end

        % Compute selected metric
        switch p_bo.ObjectiveMetric
            case "trustworthiness"
                T    = computeTrustworthinessPrecomputed(nn_high_precomputed, all_reduction, k_trust);
                loss = 1 - T;

            case "silhouette"
                sil_vals = silhouette(all_reduction(valid, :), all_labels(valid));
                loss     = -mean(sil_vals);

            case "combined"
                T        = computeTrustworthinessPrecomputed(nn_high_precomputed, all_reduction, k_trust);
                sil_vals = silhouette(all_reduction(valid, :), all_labels(valid));
                loss     = (1 - T) + (-mean(sil_vals));

            otherwise
                error('CellTypeClassifier:unknownMetric', ...
                    'Unknown ObjectiveMetric: "%s". Valid: trustworthiness, silhouette, combined.', ...
                    p_bo.ObjectiveMetric);
        end

        % Optional interneuron fraction penalty
        inh_frac = sum(all_labels(valid) == 2) / sum(valid);
        if p_bo.UseInterneuronPenalty
            loss = loss + p_bo.InterneuronWeight * (inh_frac - p_bo.InterneuronTarget)^2;
        end

        if opts.Verbose
            fprintf('  MinDist=%.3g Spread=%.3g SupervisedNDims=%d -> loss=%.3f inh=%.1f%%\n', ...
                x.MinDist, x.Spread, x.SupervisedNDims, loss, 100*inh_frac);
        end
    catch ME
        warning('CellTypeClassifier:optimizeHyperparameters', ...
            'Evaluation failed: %s', ME.message);
        loss = 1;
    end
end

% -- Run Bayesian optimization ------------------------------------------------
if opts.Verbose
    fprintf('Bayesian optimization: 3 variables (MinDist, Spread, SupervisedNDims), %d evaluations, metric=%s\n', ...
        max_evals, p_bo.ObjectiveMetric);
end

bo_result = bayesopt(@objective, vars, ...
    'MaxObjectiveEvaluations', max_evals, ...
    'IsObjectiveDeterministic', false, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'Verbose', double(opts.Verbose), ...
    'PlotFcn', {});

% -- Extract best parameters -------------------------------------------------
best = bo_result.XAtMinObjective;

results.bestParams    = struct('UMAP', struct( ...
    'MinDist',         best.MinDist, ...
    'Spread',          best.Spread, ...
    'SupervisedNDims', best.SupervisedNDims));
results.bestObjective  = bo_result.MinObjective;
results.bayesoptResult = bo_result;
results.allObjectives  = bo_result.XTrace;
results.nVars          = n_vars;

if opts.Verbose
    fprintf('\nBest topology parameters (metric=%s, %d evaluations):\n', ...
        p_bo.ObjectiveMetric, max_evals);
    fprintf('  %-25s %g\n',  'MinDist',         best.MinDist);
    fprintf('  %-25s %g\n',  'Spread',          best.Spread);
    fprintf('  %-25s %d\n',  'SupervisedNDims', best.SupervisedNDims);
    fprintf('  Best objective:         %.4f\n', results.bestObjective);
end
end

% -- Helper: trustworthiness with precomputed high-D kNN ----------------------
function T = computeTrustworthinessPrecomputed(nn_high, X_low, k)
% COMPUTETRUSTWORTHINESSPRECOMPUTED  Venna & Kaski (2006) embedding faithfulness.
%
%   T = computeTrustworthinessPrecomputed(nn_high, X_low, k)
%
%   nn_high - (n x k_hi) precomputed high-D neighbor indices (self excluded)
%   X_low   - (n x d) low-D embedding coordinates
%   k       - neighborhood size for trustworthiness evaluation
%
%   For each point, counts low-D nearest neighbors that were NOT neighbors in
%   high-D space, weighted by their high-D rank penalty.
%   T = 1 is perfect; T = 0 is worst case.
    n    = size(X_low, 1);
    k    = min(k, n - 1);
    k_hi = size(nn_high, 2);

    [nn_low, ~] = knnsearch(X_low, X_low, 'K', k + 1);
    nn_low = nn_low(:, 2:end);  % (n x k), exclude self

    penalty = 0;
    for i = 1:n
        for pos = 1:k
            j     = nn_low(i, pos);
            r_pos = find(nn_high(i, :) == j, 1);
            if isempty(r_pos)
                r_pos = k_hi + 1;
            end
            if r_pos > k
                penalty = penalty + (r_pos - k);
            end
        end
    end

    denom = n * k * (2*n - 3*k - 1);
    if denom <= 0
        T = 0;
        return
    end
    T = max(0, min(1, 1 - (2 / denom) * penalty));
end
