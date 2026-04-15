function results = optimizeUnsupervisedUMAP(ctc, opts)
% OPTIMIZEUNSUPERVISEDUMAP  Phase 1 BayOpt: unsupervised UMAP + Louvain community coherence.
%
% Optimises unsupervised UMAP hyperparameters so that drug-responsive units
% (inhibitory candidates) form a coherent Louvain community in the UMAP graph.
%
% OBJECTIVE — "community_coherence" (only option):
%   For each evaluation:
%     1. Apply ACGWeight/WaveformWeight feature group scaling.
%     2. Run unsupervised UMAP with NNeighbors, MinDist, Spread; capture umap_ev.graph.
%     3. Symmetrise umap_ev.graph (NxN fuzzy simplicial set — no kNN re-derivation).
%     4. Run community_louvain (BCT) with LouvainResolution × LouvainRestarts;
%        keep run with highest Q.
%     5. Inhibitory communities: resp_frac >= RelThresh × max_frac AND
%        resp_frac >= EnrichmentFactor × p_resp (dual criterion).
%     6. loss = -(fraction of responsive units in inhibitory communities).
%
% OPTIMIZED VARIABLES (from Parameters.BayesianOptimization.UnsupOptimizeVars):
%   NNeighbors       — unsupervised UMAP n_neighbors (integer)
%   MinDist          — UMAP min_dist (log-scaled)
%   Spread           — UMAP spread
%   ACGWeight        — ACG group L2 contribution weight
%   WaveformWeight   — Waveform group L2 contribution weight
%   LouvainResolution — Louvain gamma (continuous; higher → finer communities)
%
% After optimisation, best parameters are written to ctc.Parameters.UMAP and
% ctc.Parameters.Community.LouvainResolution. The final unsupervised UMAP is
% run once more and stored in ctc.Reduction.Unsupervised / ctc.UMAP.
%
% REQUIRES:
%   ctc.identifyResponsiveUnits() must have been run.
%   Brain Connectivity Toolbox (community_louvain) must be on the MATLAB path.
%
% USAGE:
%   ctc.identifyResponsiveUnits();
%   results = ctc.optimizeUnsupervisedUMAP();
%   ctc.Parameters = parseStructParameters(ctc.Parameters, results.bestParams);
%   ctc.generateTrainLabels();   % now uses optimised UMAP + Louvain params
%
% NAME-VALUE:
%   Verbose — print per-evaluation progress (default true)
%
% OUTPUT:
%   results.bestParams      — struct with UMAP and Community subfields
%   results.bestObjective   — best (lowest) loss = -(community_coherence)
%   results.bayesoptResult  — full BayesianOptimization object
%   results.optimizedVars   — string array of variables that were optimised
%   results.fixedVars       — string array of variables held fixed

arguments
    ctc  CellTypeClassifier
    opts.Verbose (1,1) logical = true
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before optimizeUnsupervisedUMAP()');
assert(exist('community_louvain', 'file') == 2, ...
    ['community_louvain not found on the MATLAB path. ' ...
     'Add the Brain Connectivity Toolbox (BCT) to the path first.']);

p_umap = ctc.Parameters.UMAP;
p_bo   = ctc.Parameters.BayesianOptimization;
p_comm = ctc.Parameters.Community;

% -- Build normalised feature cache (once) ------------------------------------
ctc.buildNormalizedFeatures();
nf = ctc.NormalizedFeatures;

% Global z-score on ALL unique units — matches generateTrainLabels convention.
% In the transductive setting, subset-fit statistics would bias the embedding
% toward training cultures alone. BayOpt and generateTrainLabels must use the
% same normalization so that optimized hyperparameters transfer to the main pipeline.
[X_all, ~, ~] = normalize(nf.X_pergroup);
nan_cols = any(isnan(X_all), 1);
X_all(:, nan_cols) = [];
scale_vec = max(abs(X_all), [], 1);
scale_vec(scale_vec == 0) = 1;
X_all = X_all ./ scale_vec;

% Feature selection (if enabled) — applied with current parameters
feat_names = nf.feat_names(~nan_cols);
if p_umap.FeatureSelection && size(X_all, 2) > 1
    feat_var   = var(X_all, 0, 1);
    var_thresh = prctile(feat_var, p_umap.MinVariancePercentile);
    low_var    = feat_var <= var_thresh;
    % Graph-based correlation removal (order-independent; same logic as generateTrainLabels).
    R = abs(corrcoef(X_all));
    R(logical(eye(size(R)))) = 0;
    high_corr = false(1, size(X_all, 2));
    adj = R > p_umap.MaxCorrelation;
    while any(adj(:))
        n_edges = sum(adj & ~high_corr' & ~high_corr, 2);
        n_edges(high_corr) = 0;
        [~, worst] = max(n_edges);
        high_corr(worst) = true;
        adj(worst, :) = false;
        adj(:, worst) = false;
    end
    sel_mask   = ~(low_var | high_corr);
    X_all      = X_all(:, sel_mask);
    feat_names = feat_names(sel_mask);
end

% Feature group masks
is_acg = startsWith(feat_names, "FullACG") | startsWith(feat_names, "ACG");
is_wf  = startsWith(feat_names, "Waveform") | startsWith(feat_names, "ReferenceWaveform");
n_acg  = sum(is_acg);
n_wf   = sum(is_wf);
N_all  = size(X_all, 1);

% Responsive unit flags (global unique-unit indices)
unit_ids_table = string(ctc.FeatureStore.UnitTable.UnitID);
resp_uids      = unique(unit_ids_table(ctc.ResponsiveUnitIdx));
unique_ud_ids  = string({nf.unique_ud.UnitID});
resp_unique    = ismember(unique_ud_ids, resp_uids);   % 1 × n_unique logical
n_responsive   = sum(resp_unique);
assert(n_responsive > 0, 'No responsive units found in NormalizedFeatures.');

p_resp = n_responsive / N_all;   % baseline responsive probability (for EnrichmentFactor)

if opts.Verbose
    fprintf('Phase 1 BayOpt: %d units, %d responsive (p_resp=%.3f), %d features\n', ...
        N_all, n_responsive, p_resp, size(X_all, 2));
end

% -- Fixed defaults (values used for variables not in UnsupOptimizeVars) ------
fixed_defaults = struct( ...
    'NNeighbors',        p_umap.NNeighbors, ...
    'MinDist',           p_umap.MinDist, ...
    'Spread',            p_umap.Spread, ...
    'ACGWeight',         p_umap.ACGWeight, ...
    'WaveformWeight',    p_umap.WaveformWeight, ...
    'LouvainResolution', p_comm.LouvainResolution);

% -- Adaptive LouvainResolution range -----------------------------------------
% Target: communities of ~2-5× n_responsive units. Resolution scales with
% N_all/n_responsive: fewer responsive units relative to total → need finer
% communities (higher resolution) to isolate them.
lou_range = p_bo.LouvainResolutionRange;  % hard bounds from user
if isfield(p_bo, 'AutoLouvainRange') && p_bo.AutoLouvainRange
    ratio         = N_all / max(n_responsive, 1);
    auto_res_max  = min(lou_range(2), max(lou_range(1), ratio / 10));
    lou_range_eff = [lou_range(1), auto_res_max];
    if opts.Verbose
        fprintf('  AutoLouvainRange: ratio=%.1f → LouvainResolution range [%.2f, %.2f]\n', ...
            ratio, lou_range_eff(1), lou_range_eff(2));
    end
else
    lou_range_eff = lou_range;
end

% -- Build optimizable variable array -----------------------------------------
all_valid = ["NNeighbors","MinDist","Spread","ACGWeight","WaveformWeight","LouvainResolution"];
opt_vars  = string(p_bo.UnsupOptimizeVars);
bad_vars  = opt_vars(~ismember(opt_vars, all_valid));
assert(isempty(bad_vars), ...
    'CellTypeClassifier:unknownUnsupOptVar  Unknown UnsupOptimizeVars: %s\nValid: %s', ...
    strjoin(bad_vars, ', '), strjoin(all_valid, ', '));

vars = optimizableVariable.empty;
if ismember("NNeighbors",        opt_vars)
    vars(end+1) = optimizableVariable('NNeighbors',        p_bo.NNeighborsRange,  'Type', 'integer');
end
if ismember("MinDist",           opt_vars)
    vars(end+1) = optimizableVariable('MinDist',           p_bo.MinDistRange,     'Transform', 'log');
end
if ismember("Spread",            opt_vars)
    vars(end+1) = optimizableVariable('Spread',            p_bo.SpreadRange);
end
if ismember("ACGWeight",         opt_vars)
    vars(end+1) = optimizableVariable('ACGWeight',         p_bo.ACGWeightRange);
end
if ismember("WaveformWeight",    opt_vars)
    vars(end+1) = optimizableVariable('WaveformWeight',    p_bo.WaveformWeightRange);
end
if ismember("LouvainResolution", opt_vars)
    vars(end+1) = optimizableVariable('LouvainResolution', lou_range_eff);
end

n_vars    = numel(vars);
max_evals = max(20, 10 * n_vars);
rng_seed  = ctc.Parameters.RNGSeed;

if opts.Verbose
    fprintf(['Phase 1 BayOpt: %d variables, %d evaluations\n' ...
        '  Variables: %s\n'], n_vars, max_evals, strjoin(opt_vars, ', '));
end

% -- Objective function -------------------------------------------------------
function loss = objective(x)
    v_nn  = getpval(x, 'NNeighbors',        fixed_defaults, opt_vars);
    v_md  = getpval(x, 'MinDist',           fixed_defaults, opt_vars);
    v_sp  = getpval(x, 'Spread',            fixed_defaults, opt_vars);
    v_acg = getpval(x, 'ACGWeight',         fixed_defaults, opt_vars);
    v_wf  = getpval(x, 'WaveformWeight',    fixed_defaults, opt_vars);
    v_lr  = getpval(x, 'LouvainResolution', fixed_defaults, opt_vars);

    % Apply feature group weighting
    X_ev = X_all;
    if n_acg > 0; X_ev(:, is_acg) = X_ev(:, is_acg) * sqrt(v_acg / n_acg); end
    if n_wf  > 0; X_ev(:, is_wf)  = X_ev(:, is_wf)  * sqrt(v_wf  / n_wf);  end

    % Run unsupervised UMAP — capture model object for its graph
    rng(rng_seed, 'twister');
    k_nn = min(v_nn, N_all - 1);
    try
        [~, umap_ev, ~, ~] = run_umap(X_ev, ...
            'n_components', p_umap.NDims, ...
            'n_neighbors',  k_nn, ...
            'min_dist',     v_md, ...
            'spread',       v_sp, ...
            'method', 'Java', 'verbose', 'none');
    catch ME
        warning('CellTypeClassifier:optimizeUnsupervisedUMAP', ...
            'UMAP failed: %s', ME.message);
        loss = 1;
        return
    end

    if isempty(umap_ev) || isempty(umap_ev.head)
        loss = 1; return
    end

    % Reconstruct weighted graph from COO edge list stored by simplicial_set_embedding.
    % umap.graph is not reliably populated after a Java SGD run; head/tail/weights
    % are always set and carry the same fuzzy membership strengths (pruned of very
    % weak edges < max_weight/n_epochs, but fully weighted otherwise).
    G = sparse(umap_ev.head, umap_ev.tail, umap_ev.weights, N_all, N_all);
    G = (G + G') / 2;

    % Run community_louvain with multiple restarts; keep highest Q
    W_lou = full(G);   % BCT requires full matrix
    best_Q  = -Inf;
    best_M  = ones(N_all, 1);
    n_restarts = p_comm.LouvainRestarts;
    for r = 1:n_restarts
        rng(rng_seed + r, 'twister');
        try
            [M_r, Q_r] = community_louvain(W_lou, v_lr);
            if Q_r > best_Q
                best_Q = Q_r;
                best_M = M_r;
            end
        catch
            % individual restart failed — continue
        end
    end

    if best_Q == -Inf
        loss = 1; return   % all restarts failed
    end

    % Identify inhibitory communities (high responsive fraction)
    n_comm = max(best_M);
    comm_resp_frac = zeros(1, n_comm);
    for c = 1:n_comm
        c_mask = best_M == c;
        if any(c_mask)
            comm_resp_frac(c) = sum(resp_unique(c_mask)) / sum(c_mask);
        end
    end
    max_frac = max(comm_resp_frac);
    if max_frac == 0
        loss = 1; return
    end
    inh_by_rel = comm_resp_frac >= max_frac * p_comm.InhibitoryCommunityRelThresh;
    inh_by_abs = comm_resp_frac >= p_comm.EnrichmentFactor * p_resp;
    inh_mask   = inh_by_rel & inh_by_abs;

    % Community coherence: fraction of responsive units in inhibitory communities
    resp_in_inh = sum(resp_unique & inh_mask(best_M)');
    coherence   = resp_in_inh / n_responsive;
    loss        = -coherence;

    if opts.Verbose
        fprintf(['  NNeigh=%d MinDist=%.3g Spread=%.3g ACG=%.2f WF=%.2f ' ...
            'LouRes=%.2f → coherence=%.3f Q=%.3f n_comm=%d\n'], ...
            v_nn, v_md, v_sp, v_acg, v_wf, v_lr, coherence, best_Q, n_comm);
    end
end

% -- Run Bayesian optimisation ------------------------------------------------
bo_result = bayesopt(@objective, vars, ...
    'MaxObjectiveEvaluations', max_evals, ...
    'IsObjectiveDeterministic', false, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'Verbose', double(opts.Verbose), ...
    'PlotFcn', {});

% -- Extract best parameters --------------------------------------------------
best      = bo_result.XAtMinObjective;
best_full = fixed_defaults;
for vn = opt_vars
    best_full.(vn) = best.(vn);
end

results.bestParams = struct( ...
    'UMAP', struct( ...
        'NNeighbors',   best_full.NNeighbors, ...
        'MinDist',      best_full.MinDist, ...
        'Spread',       best_full.Spread, ...
        'ACGWeight',    best_full.ACGWeight, ...
        'WaveformWeight', best_full.WaveformWeight), ...
    'Community', struct( ...
        'LouvainResolution', best_full.LouvainResolution));
results.bestObjective  = bo_result.MinObjective;
results.bayesoptResult = bo_result;
results.optimizedVars  = opt_vars;
results.fixedVars      = all_valid(~ismember(all_valid, opt_vars));

if opts.Verbose
    fprintf('\nBest Phase 1 parameters (%d evaluations):\n', max_evals);
    all_var_names = ["NNeighbors","MinDist","Spread","ACGWeight","WaveformWeight","LouvainResolution"];
    for vn = all_var_names
        tag = '';
        if ~ismember(vn, opt_vars); tag = ' (fixed)'; end
        fprintf('  %-25s %g%s\n', vn, best_full.(vn), tag);
    end
    fprintf('  Best coherence:         %.4f\n', -results.bestObjective);
end

% -- Apply best params and run final unsupervised UMAP ------------------------
ctc.Parameters.UMAP.NNeighbors    = best_full.NNeighbors;
ctc.Parameters.UMAP.MinDist       = best_full.MinDist;
ctc.Parameters.UMAP.Spread        = best_full.Spread;
ctc.Parameters.UMAP.ACGWeight     = best_full.ACGWeight;
ctc.Parameters.UMAP.WaveformWeight= best_full.WaveformWeight;
ctc.Parameters.Community.LouvainResolution = best_full.LouvainResolution;

if opts.Verbose
    fprintf('Running final unsupervised UMAP with best parameters...\n');
end
% Feature weighting for final run
X_final = X_all;
if n_acg > 0; X_final(:, is_acg) = X_final(:, is_acg) * sqrt(best_full.ACGWeight  / n_acg); end
if n_wf  > 0; X_final(:, is_wf)  = X_final(:, is_wf)  * sqrt(best_full.WaveformWeight / n_wf);  end

rng(rng_seed, 'twister');
[reduction_final, umap_model, ~, ~] = run_umap(X_final, ...
    'n_components', p_umap.NDims, ...
    'n_neighbors',  min(best_full.NNeighbors, N_all - 1), ...
    'min_dist',     best_full.MinDist, ...
    'spread',       best_full.Spread, ...
    'method', 'Java', 'verbose', 'none');
ctc.UMAP = umap_model;
if isempty(ctc.Reduction); ctc.Reduction = struct(); end
ctc.Reduction.Unsupervised = reduction_final;
fprintf('Phase 1 complete. Final UMAP stored in ctc.Reduction.Unsupervised.\n');

if ctc.Parameters.Diagnostics.Enable
    ctc.diagnosticOptimization(results);
end
end

% -- Helper: resolve parameter value (optimised table or fixed default) -------
function val = getpval(x, name, fixed_defaults, opt_vars)
    if ismember(name, opt_vars)
        val = x.(name);
    else
        val = fixed_defaults.(name);
    end
end
