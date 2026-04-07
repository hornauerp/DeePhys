function results = optimizeHyperparameters(ctc, opts)
% OPTIMIZEHYPERPARAMETERS  Bayesian optimization of UMAP parameters.
%
% Optimizes 6 UMAP parameters to produce well-structured supervised embeddings.
%
% OBJECTIVE METRIC (set via Parameters.BayesianOptimization.ObjectiveMetric):
%   "weighted_silhouette" (default) — evaluates cluster separation on the TEST
%       embedding only. Asymmetric weighting: -(0.7*inh_sil + 0.3*exc_sil).
%       Tested on test units to avoid conflating parameter quality with the
%       strength of the supervised signal in the train embedding.
%   "trustworthiness" — 1 - T, where T (Venna & Kaski 2006) measures how
%       faithfully the low-D embedding preserves high-D neighborhoods.
%   "silhouette"       — negative mean silhouette on supervised embedding.
%   "combined"         — (1-T) + (-mean_sil).
%
% OPTIMIZED VARIABLES (6):
%   MinDist              — UMAP min_dist (log-scaled over MinDistRange)
%   Spread               — UMAP spread (over SpreadRange)
%   SupervisedNDims      — supervised UMAP output dimensions (integer)
%   SupervisedNNeighbors — n_neighbors for supervised UMAP (integer, ceiling applied)
%   TargetWeight         — label vs data topology balance
%   ACGWeight            — ACG group L2 contribution weight
%   WaveformWeight       — Waveform group L2 contribution weight
%
% Fixed: 60 Bayesian evaluations (10 per variable).
%
% Efficiency: training labels and normalised feature matrix are generated once.
% ACGWeight/WaveformWeight are applied per-evaluation before supervisedUMAP.
% High-D kNN is precomputed once for trustworthiness (if that metric is used).
%
% Requires identifyResponsiveUnits() to have been run first.
%
% USAGE:
%   ctc.identifyResponsiveUnits();
%   results = ctc.optimizeHyperparameters();
%   ctc.Parameters = parseStructParameters(ctc.Parameters, results.bestParams);
%   ctc.generateTrainLabels();
%   ctc.classifyUnits();
%
% NAME-VALUE:
%   Verbose - Print progress during optimization (default true)
%
% OUTPUTS:
%   results - struct with fields:
%     .bestParams     - struct of optimal parameter values (UMAP subfield)
%     .bestObjective  - best (lowest) objective value achieved
%     .bayesoptResult - full BayesianOptimization object
%     .allObjectives  - table of all evaluations
%     .nVars          - number of variables optimized (7)

arguments
    ctc CellTypeClassifier
    opts.Verbose (1,1) logical = true
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before optimizeHyperparameters()');

p_bo = ctc.Parameters.BayesianOptimization;

% -- Define 6 optimizable variables ------------------------------------------
vars = [ ...
    optimizableVariable('MinDist',              p_bo.MinDistRange,              'Transform', 'log'), ...
    optimizableVariable('Spread',               p_bo.SpreadRange), ...
    optimizableVariable('SupervisedNDims',      p_bo.SupervisedNDimsRange,      'Type', 'integer'), ...
    optimizableVariable('SupervisedNNeighbors', p_bo.SupervisedNNeighborsRange, 'Type', 'integer'), ...
    optimizableVariable('TargetWeight',         p_bo.TargetWeightRange), ...
    optimizableVariable('ACGWeight',            p_bo.ACGWeightRange), ...
    optimizableVariable('WaveformWeight',       p_bo.WaveformWeightRange)];
n_vars    = 7;
max_evals = 60;

% -- Generate training labels once (fixed across all iterations) --------------
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

% -- Precompute base feature matrix (before per-eval weighting) ---------------
X_base = normalize(fixed_nf.X_pergroup, 'center', fixed_np.mu_global, 'scale', fixed_np.sigma_global);
X_base(:, fixed_np.nan_cols) = [];
X_base = X_base ./ fixed_np.scale;
if isfield(fixed_np, 'feature_selection_mask') && ~all(fixed_np.feature_selection_mask)
    X_base = X_base(:, fixed_np.feature_selection_mask);
end

% Feature group masks for weighting
feat_names_trimmed = fixed_np.feat_names_trimmed;
is_acg = startsWith(feat_names_trimmed, "FullACG") | startsWith(feat_names_trimmed, "ACG");
is_wf  = startsWith(feat_names_trimmed, "Waveform") | startsWith(feat_names_trimmed, "ReferenceWaveform");
n_acg  = sum(is_acg);
n_wf   = sum(is_wf);

% Minority class ceiling constant (same formula as classifyUnits)
n_inh_train = sum(fixed_labels.sorted_y_train == 2);
minority_ceiling = max(1, floor(n_inh_train / 5));

train_mask = logical(fixed_labels.umap_train_idx);

% -- Precompute high-D kNN for trustworthiness (if needed) --------------------
% Only computed when ObjectiveMetric involves trustworthiness — expensive step.
nn_high_precomputed = [];
k_trust = 15;
if ismember(p_bo.ObjectiveMetric, ["trustworthiness", "combined"])
    X_all_ordered = [X_base(train_mask, :); X_base(~train_mask, :)];
    k_hi    = min(5 * k_trust, size(X_all_ordered, 1) - 1);
    [nn_high_precomputed, ~] = knnsearch(X_all_ordered, X_all_ordered, 'K', k_hi + 1);
    nn_high_precomputed = nn_high_precomputed(:, 2:end);
    if opts.Verbose
        fprintf('Precomputed high-D kNN (k_hi=%d) for trustworthiness.\n', k_hi);
    end
end

if opts.Verbose
    fprintf('Precomputed: %d unique units (%d train), %d features\n', ...
        n_unique, sum(train_mask), size(X_base, 2));
end

% -- Objective function -------------------------------------------------------
function loss = objective(x)
    % Apply per-evaluation feature group weighting
    X_eval = X_base;
    if n_acg > 0
        X_eval(:, is_acg) = X_eval(:, is_acg) * sqrt(x.ACGWeight / n_acg);
    end
    if n_wf > 0
        X_eval(:, is_wf) = X_eval(:, is_wf) * sqrt(x.WaveformWeight / n_wf);
    end

    % Apply minority ceiling to candidate SupervisedNNeighbors
    sup_nn = min(x.SupervisedNNeighbors, minority_ceiling);

    temp_params = ctc.Parameters;
    temp_params.UMAP.MinDist              = x.MinDist;
    temp_params.UMAP.Spread               = x.Spread;
    temp_params.UMAP.SupervisedNDims      = x.SupervisedNDims;
    temp_params.UMAP.SupervisedNNeighbors = sup_nn;
    temp_params.UMAP.TargetWeight         = x.TargetWeight;
    temp_params.UMAP.ACGWeight            = x.ACGWeight;
    temp_params.UMAP.WaveformWeight       = x.WaveformWeight;
    temp_params.Diagnostics.Enable        = false;  % suppress diagnostics inside BayOpt

    temp_ctc = CellTypeClassifier(ctc.FeatureStore, ctc.UnitDataArray, temp_params);
    temp_ctc.TrainLabels         = fixed_labels;
    temp_ctc.NormalizationParams = fixed_np;
    temp_ctc.NormalizedFeatures  = fixed_nf;
    temp_ctc.Reduction           = struct();

    % Override the normalized feature cache so classifyUnits uses our pre-weighted X
    % We do this by injecting a modified NormalizedFeatures where X_pergroup = X_eval
    % and blanking the normalization params so classifyUnits applies identity transform
    nf_override            = fixed_nf;
    nf_override.X_pergroup = X_eval;
    temp_ctc.NormalizedFeatures = nf_override;

    % Override NormalizationParams to identity (weighting already applied above)
    np_identity                         = fixed_np;
    np_identity.mu_global               = zeros(1, size(X_eval, 2));
    np_identity.sigma_global            = ones(1, size(X_eval, 2));
    np_identity.nan_cols                = false(1, size(X_eval, 2));
    np_identity.scale                   = ones(1, size(X_eval, 2));
    np_identity.feature_selection_mask  = true(1, size(X_eval, 2));
    % feat_names_trimmed preserved so weighting block in classifyUnits is a no-op
    % (ACGWeight=1/n_acg * sqrt(n_acg) = 1 at the already-scaled level)
    temp_ctc.NormalizationParams = np_identity;

    try
        temp_ctc.classifyUnits();

        all_reduction = [temp_ctc.Reduction.Train; temp_ctc.Reduction.Test];
        unique_labels = temp_ctc.UnitLabels(fixed_nf.unique_to_rep);
        train_labels_vec = unique_labels(train_mask)';
        test_labels_vec  = unique_labels(~train_mask)';
        all_labels       = [train_labels_vec; test_labels_vec];
        valid            = ~isnan(all_labels);

        if sum(valid) < 10 || numel(unique(all_labels(valid))) < 2
            loss = 1;
            return
        end

        switch p_bo.ObjectiveMetric
            case "weighted_silhouette"
                % Evaluate on test embedding only
                valid_test = ~isnan(test_labels_vec);
                if sum(valid_test) < 10 || numel(unique(test_labels_vec(valid_test))) < 2
                    loss = 1; return
                end
                test_red = temp_ctc.Reduction.Test;
                sil_vals = silhouette(test_red(valid_test, :), test_labels_vec(valid_test));
                inh_mask_t = test_labels_vec(valid_test) == 2;
                exc_mask_t = test_labels_vec(valid_test) == 1;
                inh_sil = 0;
                exc_sil = 0;
                if any(inh_mask_t); inh_sil = mean(sil_vals(inh_mask_t)); end
                if any(exc_mask_t); exc_sil = mean(sil_vals(exc_mask_t)); end
                loss = -(0.7 * inh_sil + 0.3 * exc_sil);

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
                    'Unknown ObjectiveMetric: "%s". Valid: weighted_silhouette, trustworthiness, silhouette, combined.', ...
                    p_bo.ObjectiveMetric);
        end

        % Optional interneuron fraction penalty
        inh_frac = sum(all_labels(valid) == 2) / sum(valid);
        if p_bo.UseInterneuronPenalty
            loss = loss + p_bo.InterneuronWeight * (inh_frac - p_bo.InterneuronTarget)^2;
        end

        if opts.Verbose
            fprintf(['  MinDist=%.3g Spread=%.3g NDims=%d NNeigh=%d(->%d) ' ...
                'TW=%.3f ACG=%.2f WF=%.2f -> loss=%.3f inh=%.1f%%\n'], ...
                x.MinDist, x.Spread, x.SupervisedNDims, x.SupervisedNNeighbors, sup_nn, ...
                x.TargetWeight, x.ACGWeight, x.WaveformWeight, loss, 100*inh_frac);
        end
    catch ME
        warning('CellTypeClassifier:optimizeHyperparameters', ...
            'Evaluation failed: %s', ME.message);
        loss = 1;
    end
end

% -- Run Bayesian optimization ------------------------------------------------
if opts.Verbose
    fprintf(['Bayesian optimization: %d variables, %d evaluations, metric=%s\n' ...
        '  Variables: MinDist, Spread, SupervisedNDims, SupervisedNNeighbors, ' ...
        'TargetWeight, ACGWeight, WaveformWeight\n'], n_vars, max_evals, p_bo.ObjectiveMetric);
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
    'MinDist',              best.MinDist, ...
    'Spread',               best.Spread, ...
    'SupervisedNDims',      best.SupervisedNDims, ...
    'SupervisedNNeighbors', best.SupervisedNNeighbors, ...
    'TargetWeight',         best.TargetWeight, ...
    'ACGWeight',            best.ACGWeight, ...
    'WaveformWeight',       best.WaveformWeight));
results.bestObjective  = bo_result.MinObjective;
results.bayesoptResult = bo_result;
results.allObjectives  = bo_result.XTrace;
results.nVars          = n_vars;

if opts.Verbose
    fprintf('\nBest parameters (metric=%s, %d evaluations):\n', p_bo.ObjectiveMetric, max_evals);
    fprintf('  %-25s %g\n',  'MinDist',              best.MinDist);
    fprintf('  %-25s %g\n',  'Spread',               best.Spread);
    fprintf('  %-25s %d\n',  'SupervisedNDims',      best.SupervisedNDims);
    fprintf('  %-25s %d\n',  'SupervisedNNeighbors', best.SupervisedNNeighbors);
    fprintf('  %-25s %.3f\n','TargetWeight',          best.TargetWeight);
    fprintf('  %-25s %.3f\n','ACGWeight',             best.ACGWeight);
    fprintf('  %-25s %.3f\n','WaveformWeight',        best.WaveformWeight);
    fprintf('  Best objective:         %.4f\n', results.bestObjective);
end

% -- Optional diagnostic ------------------------------------------------------
if ctc.Parameters.Diagnostics.Enable
    ctc.diagnosticOptimization(results);
end
end

% -- Helper: trustworthiness with precomputed high-D kNN ----------------------
function T = computeTrustworthinessPrecomputed(nn_high, X_low, k)
% COMPUTETRUSTWORTHINESSPRECOMPUTED  Venna & Kaski (2006) embedding faithfulness.
    n    = size(X_low, 1);
    k    = min(k, n - 1);
    k_hi = size(nn_high, 2);

    [nn_low, ~] = knnsearch(X_low, X_low, 'K', k + 1);
    nn_low = nn_low(:, 2:end);

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
