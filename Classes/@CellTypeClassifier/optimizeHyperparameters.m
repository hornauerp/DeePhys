function results = optimizeHyperparameters(ctc, opts)
% OPTIMIZEHYPERPARAMETERS  Bayesian optimization of UMAP parameters.
%
% Optimizes 6 UMAP parameters to produce well-structured supervised embeddings.
%
% OBJECTIVE METRIC (set via Parameters.BayesianOptimization.ObjectiveMetric):
%   "edge_homophily" (default) — fraction of k-NN in the supervised embedding
%       that share the same predicted label. Directly aligns Phase 2 BayOpt with
%       the graph classifier: a higher-homophily embedding makes graph label
%       propagation more reliable. Evaluated on test units only (pseudo-labelled
%       by kNN in feature space) to avoid trivial train-node anchoring.
%   "weighted_silhouette" — evaluates cluster separation on the TEST
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

warning('CellTypeClassifier:deprecated', ...
    ['optimizeHyperparameters (Phase 2) is deprecated. ' ...
     'Supervised UMAP parameters are no longer used in classification — ' ...
     'classifyUnits() propagates labels on the unsupervised ctc.UMAP.graph. ' ...
     'Use optimizeUnsupervisedUMAP() for Phase 1 hyperparameter optimization.']);

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before optimizeHyperparameters()');

p_bo = ctc.Parameters.BayesianOptimization;

% -- Fixed defaults: values used for variables not being optimized ------------
% Sourced from current ctc.Parameters so the user can pre-set them before calling.
fixed_defaults = struct( ...
    'MinDist',              ctc.Parameters.UMAP.MinDist, ...
    'Spread',               ctc.Parameters.UMAP.Spread, ...
    'SupervisedNDims',      ctc.Parameters.UMAP.SupervisedNDims, ...
    'SupervisedNNeighbors', ctc.Parameters.UMAP.SupervisedNNeighbors, ...
    'TargetWeight',         ctc.Parameters.UMAP.TargetWeight, ...
    'ACGWeight',            ctc.Parameters.UMAP.ACGWeight, ...
    'WaveformWeight',       ctc.Parameters.UMAP.WaveformWeight);

% -- Build optimizable variable array from OptimizeVars -----------------------
all_valid = ["MinDist", "Spread", "SupervisedNDims", "SupervisedNNeighbors", ...
             "TargetWeight", "ACGWeight", "WaveformWeight"];
opt_vars  = string(p_bo.OptimizeVars);
bad_vars  = opt_vars(~ismember(opt_vars, all_valid));
assert(isempty(bad_vars), ...
    'CellTypeClassifier:unknownOptVar  Unknown OptimizeVars: %s\nValid: %s', ...
    strjoin(bad_vars, ', '), strjoin(all_valid, ', '));

vars = optimizableVariable.empty;
if ismember("MinDist",              opt_vars)
    vars(end+1) = optimizableVariable('MinDist',              p_bo.MinDistRange,              'Transform', 'log');
end
if ismember("Spread",               opt_vars)
    vars(end+1) = optimizableVariable('Spread',               p_bo.SpreadRange);
end
if ismember("SupervisedNDims",      opt_vars)
    vars(end+1) = optimizableVariable('SupervisedNDims',      p_bo.SupervisedNDimsRange,      'Type', 'integer');
end
if ismember("SupervisedNNeighbors", opt_vars)
    vars(end+1) = optimizableVariable('SupervisedNNeighbors', p_bo.SupervisedNNeighborsRange, 'Type', 'integer');
end
if ismember("TargetWeight",         opt_vars)
    vars(end+1) = optimizableVariable('TargetWeight',         p_bo.TargetWeightRange);
end
if ismember("ACGWeight",            opt_vars)
    vars(end+1) = optimizableVariable('ACGWeight',            p_bo.ACGWeightRange);
end
if ismember("WaveformWeight",       opt_vars)
    vars(end+1) = optimizableVariable('WaveformWeight',       p_bo.WaveformWeightRange);
end

n_vars    = numel(vars);
max_evals = max(20, 10 * n_vars);

% -- Generate training labels once (fixed across all iterations) --------------
if opts.Verbose
    fprintf('Generating fixed training labels for optimization...\n');
end
ref_ctc = CellTypeClassifier(ctc.FeatureStore, ctc.UnitDataArray, ctc.Parameters);
ref_ctc.ResponsiveUnitIdx       = ctc.ResponsiveUnitIdx;
ref_ctc.ResponsiveUnitDirection = ctc.ResponsiveUnitDirection;
ref_ctc.ResponsiveStrength      = ctc.ResponsiveStrength;
ref_ctc.CounterexampleUnitIdx   = ctc.CounterexampleUnitIdx;
ref_ctc.OriginalUnitDataArray   = ctc.OriginalUnitDataArray;
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

% -- Precompute train/test split + kNN settings (avoids per-eval overhead) ----
X_train_base = X_base(train_mask, :);
X_test_base  = X_base(~train_mask, :);
y_train      = fixed_labels.sorted_y_train;
conf_k       = ctc.Parameters.UMAP.ConfidenceK;
if ctc.Parameters.UMAP.AutoConfidenceK
    conf_k = max(5, round(sqrt(sum(train_mask))));
end
tdir = tempdir;   % ephemeral templates — faster than user TemplateDir
rng_seed = ctc.Parameters.RNGSeed;

if opts.Verbose
    fprintf('Precomputed: %d unique units (%d train), %d features\n', ...
        n_unique, sum(train_mask), size(X_base, 2));
end

% -- Objective function (inlined UMAP + kNN; no CellTypeClassifier per eval) --
function loss = objective(x)
    % Resolve each parameter: from BayOpt table if being optimised, else fixed default
    v_min_dist  = getpval(x, 'MinDist',              fixed_defaults, opt_vars);
    v_spread    = getpval(x, 'Spread',               fixed_defaults, opt_vars);
    v_ndims     = getpval(x, 'SupervisedNDims',      fixed_defaults, opt_vars);
    v_nneigh    = getpval(x, 'SupervisedNNeighbors', fixed_defaults, opt_vars);
    v_tw        = getpval(x, 'TargetWeight',         fixed_defaults, opt_vars);
    v_acg       = getpval(x, 'ACGWeight',            fixed_defaults, opt_vars);
    v_wf        = getpval(x, 'WaveformWeight',       fixed_defaults, opt_vars);

    % Apply per-evaluation feature group weighting
    X_tr = X_train_base;
    X_te = X_test_base;
    if n_acg > 0
        X_tr(:, is_acg) = X_tr(:, is_acg) * sqrt(v_acg / n_acg);
        X_te(:, is_acg) = X_te(:, is_acg) * sqrt(v_acg / n_acg);
    end
    if n_wf > 0
        X_tr(:, is_wf) = X_tr(:, is_wf) * sqrt(v_wf / n_wf);
        X_te(:, is_wf) = X_te(:, is_wf) * sqrt(v_wf / n_wf);
    end

    sup_nn = min(v_nneigh, minority_ceiling);

    % -- Run supervised UMAP directly (training + test projection) --------
    rng(rng_seed, 'twister');
    uid           = char(java.util.UUID.randomUUID);
    template_file = fullfile(tdir, ['ctc_bo_template_', uid, '.mat']);
    clean_up      = onCleanup(@() deleteFile(template_file));

    try
        [train_red, ~, ~, ~] = run_umap([X_tr, y_train'], ...
            'label_column',       'end', ...
            'n_components',       v_ndims, ...
            'n_neighbors',        sup_nn, ...
            'min_dist',           v_min_dist, ...
            'spread',             v_spread, ...
            'target_weight',      v_tw, ...
            ...
            'method',             'java', ...
            'verbose',            'none', ...
            'save_template_file', template_file);

        [test_red, ~, ~, ~] = run_umap(X_te, ...
            'method',            'java', ...
            'verbose',           'none', ...
            'match_supervisors', 1, ...
            'template_file',     template_file);
    catch ME
        warning('CellTypeClassifier:optimizeHyperparameters', ...
            'UMAP failed: %s', ME.message);
        loss = 1;
        return
    end

    % -- kNN classification in feature space (not embedding) -------------
    % test_labels are stable across BayOpt evals; only the embedding-based
    % loss metric (silhouette/trustworthiness on test_red) varies.
    k_actual = min(conf_k, size(X_tr, 1));
    [nn_idx, nn_dists] = knnsearch(X_tr, X_te, 'K', k_actual);
    n_test      = size(test_red, 1);
    test_labels = nan(n_test, 1);
    for i = 1:n_test
        nb    = y_train(nn_idx(i, :));
        d     = nn_dists(i, :);
        eps_d = max(d) * 1e-6;
        if eps_d == 0; eps_d = 1e-10; end
        w = 1 ./ (d + eps_d);
        if sum(w(nb == 1)) >= sum(w(nb == 2))
            test_labels(i) = 1;
        else
            test_labels(i) = 2;
        end
    end

    all_labels = [y_train(:); test_labels];
    valid      = ~isnan(all_labels);

    if sum(valid) < 10 || numel(unique(all_labels(valid))) < 2
        loss = 1;
        return
    end

    % -- Objective metric -------------------------------------------------
    switch p_bo.ObjectiveMetric
        case "edge_homophily"
            % Fraction of k-NN in supervised embedding sharing the same
            % predicted label. Evaluated on test nodes only (train labels
            % are fixed anchors and would trivially inflate homophily).
            all_red     = [train_red; test_red];
            N_all_h     = size(all_red, 1);
            N_tr_h      = numel(y_train);
            k_hom       = min(conf_k, N_all_h - 1);
            [nn_emb, ~] = knnsearch(all_red, all_red, 'K', k_hom + 1);
            nn_emb      = nn_emb(:, 2:end);   % (N_all_h x k_hom)
            all_lbl     = [y_train(:); test_labels(:)];
            nn_lbl      = reshape(all_lbl(nn_emb(:)), N_all_h, k_hom);
            hom_all     = mean(nn_lbl == all_lbl, 2);   % (N_all_h x 1)
            test_hom    = hom_all(N_tr_h + 1:end);
            valid_th    = ~isnan(test_labels);
            if sum(valid_th) < 5
                loss = 1; return
            end
            loss = -mean(test_hom(valid_th));

        case "weighted_silhouette"
            valid_test = ~isnan(test_labels);
            if sum(valid_test) < 10 || numel(unique(test_labels(valid_test))) < 2
                loss = 1; return
            end
            sil_vals   = silhouette(test_red(valid_test, :), test_labels(valid_test));
            inh_mask_t = test_labels(valid_test) == 2;
            exc_mask_t = test_labels(valid_test) == 1;
            inh_sil = 0; exc_sil = 0;
            if any(inh_mask_t); inh_sil = mean(sil_vals(inh_mask_t)); end
            if any(exc_mask_t); exc_sil = mean(sil_vals(exc_mask_t)); end
            loss = -(0.7 * inh_sil + 0.3 * exc_sil);

        case "trustworthiness"
            all_red = [train_red; test_red];
            T    = computeTrustworthinessPrecomputed(nn_high_precomputed, all_red, k_trust);
            loss = 1 - T;

        case "silhouette"
            all_red  = [train_red; test_red];
            sil_vals = silhouette(all_red(valid, :), all_labels(valid));
            loss     = -mean(sil_vals);

        case "combined"
            all_red  = [train_red; test_red];
            T        = computeTrustworthinessPrecomputed(nn_high_precomputed, all_red, k_trust);
            sil_vals = silhouette(all_red(valid, :), all_labels(valid));
            loss     = (1 - T) + (-mean(sil_vals));

        otherwise
            error('CellTypeClassifier:unknownMetric', ...
                'Unknown ObjectiveMetric: "%s". Valid: edge_homophily, weighted_silhouette, trustworthiness, silhouette, combined.', ...
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
            v_min_dist, v_spread, v_ndims, v_nneigh, sup_nn, ...
            v_tw, v_acg, v_wf, loss, 100*inh_frac);
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
% Merge optimised values with fixed defaults so bestParams is always complete.
best = bo_result.XAtMinObjective;
best_full = fixed_defaults;
for vn = opt_vars
    best_full.(vn) = best.(vn);
end

results.bestParams    = struct('UMAP', struct( ...
    'MinDist',              best_full.MinDist, ...
    'Spread',               best_full.Spread, ...
    'SupervisedNDims',      best_full.SupervisedNDims, ...
    'SupervisedNNeighbors', best_full.SupervisedNNeighbors, ...
    'TargetWeight',         best_full.TargetWeight, ...
    'ACGWeight',            best_full.ACGWeight, ...
    'WaveformWeight',       best_full.WaveformWeight));
results.bestObjective  = bo_result.MinObjective;
results.bayesoptResult = bo_result;
results.allObjectives  = bo_result.XTrace;
results.nVars          = n_vars;
results.optimizedVars  = opt_vars;
results.fixedVars      = all_valid(~ismember(all_valid, opt_vars));

if opts.Verbose
    fprintf('\nBest parameters (metric=%s, %d evaluations):\n', p_bo.ObjectiveMetric, max_evals);
    all_var_names = ["MinDist","Spread","SupervisedNDims","SupervisedNNeighbors", ...
                     "TargetWeight","ACGWeight","WaveformWeight"];
    for vn = all_var_names
        tag = '';
        if ~ismember(vn, opt_vars); tag = ' (fixed)'; end
        fprintf('  %-25s %g%s\n', vn, best_full.(vn), tag);
    end
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

% -- Helper: resolve parameter value (optimised table or fixed default) -------
function val = getpval(x, name, fixed_defaults, opt_vars)
% Returns x.(name) if 'name' is in opt_vars, otherwise fixed_defaults.(name).
    if ismember(name, opt_vars)
        val = x.(name);
    else
        val = fixed_defaults.(name);
    end
end

% -- Helper: safe file deletion -----------------------------------------------
function deleteFile(f)
if isfile(f); delete(f); end
end
