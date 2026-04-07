function results = optimizeHyperparameters(ctc, opts)
% OPTIMIZEHYPERPARAMETERS  Bayesian optimization of UMAP topology parameters.
%
% Optimizes the UMAP embedding geometry parameters (MinDist, Spread,
% NNeighbors) to produce a well-structured supervised embedding. These are
% the parameters that directly control embedding shape — fixing problems
% like degenerate ellipses or collapsed clusters.
%
% Only topology parameters are optimized. Parameters that affect label
% assignment or class balance (TargetWeight, CounterexampleRatio) are
% excluded because the objective metric (QF dissimilarity / silhouette)
% can be trivially gamed by these — e.g. high TargetWeight always produces
% tighter clusters regardless of label quality.
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
%   "qf_dissimilarity" (default) — 1 - QF overlap from UMAP toolbox.
%       Measures topology-label agreement: did UMAP find real structure
%       that aligns with the class labels? Not gameable by topology params.
%   "silhouette"       — negative mean silhouette on supervised embedding.
%   "combined"         — negative silhouette + (1 - QF overlap/100).
%
% OPTIMIZED VARIABLES (3-5, depending on Auto* flags):
%   MinDist              — UMAP min_dist (always included, log-scaled)
%   Spread               — UMAP spread (always included)
%   NNeighbors           — unsupervised UMAP n_neighbors (excluded if AutoNNeighbors)
%   SupervisedNNeighbors — supervised UMAP n_neighbors (excluded if AutoNNeighbors)
%   NDims                — unsupervised UMAP output dimensions (excluded if AutoSupervisedNDims)
%
% OUTPUTS:
%   results - struct with fields:
%     .bestParams     - struct of optimal parameter values
%     .bestObjective  - best (lowest) objective value achieved
%     .bayesoptResult - full BayesianOptimization object
%     .allObjectives  - table of all evaluations
%     .nVars          - number of variables optimized

arguments
    ctc CellTypeClassifier
    opts.Verbose (1,1) logical = true
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before optimizeHyperparameters()');

p_bo   = ctc.Parameters.BayesianOptimization;
p_umap = ctc.Parameters.UMAP;

% -- Define optimizable variables: topology parameters only ------------------
vars = optimizableVariable.empty;
var_names = {};

% MinDist and Spread always included — primary geometry controls
vars(end+1) = optimizableVariable('MinDist', p_bo.MinDistRange, 'Transform', 'log');
var_names{end+1} = 'MinDist';
vars(end+1) = optimizableVariable('Spread', p_bo.SpreadRange);
var_names{end+1} = 'Spread';

% NNeighbors: excluded if AutoNNeighbors handles it
if ~p_umap.AutoNNeighbors
    vars(end+1) = optimizableVariable('NNeighbors', p_bo.NNeighborsRange, 'Type', 'integer');
    var_names{end+1} = 'NNeighbors';
    vars(end+1) = optimizableVariable('SupervisedNNeighbors', p_bo.SupervisedNNeighborsRange, 'Type', 'integer');
    var_names{end+1} = 'SupervisedNNeighbors';
end

% NDims: excluded if AutoSupervisedNDims handles it via TWO-NN
if ~p_umap.AutoSupervisedNDims
    vars(end+1) = optimizableVariable('NDims', p_bo.NDimsRange, 'Type', 'integer');
    var_names{end+1} = 'NDims';
end

n_vars    = numel(vars);
max_evals = max(15, 10 * n_vars);

% -- Objective function -------------------------------------------------------
function loss = objective(x)
    temp_params = ctc.Parameters;

    % Apply candidate topology values
    temp_params.UMAP.MinDist = x.MinDist;
    temp_params.UMAP.Spread  = x.Spread;
    if ismember('NNeighbors',          var_names); temp_params.UMAP.NNeighbors          = x.NNeighbors; end
    if ismember('SupervisedNNeighbors',var_names); temp_params.UMAP.SupervisedNNeighbors = x.SupervisedNNeighbors; end
    if ismember('NDims',               var_names); temp_params.UMAP.NDims               = x.NDims; end

    temp_ctc = CellTypeClassifier(ctc.FeatureStore, ctc.UnitDataArray, temp_params);
    temp_ctc.ResponsiveUnitIdx       = ctc.ResponsiveUnitIdx;
    temp_ctc.ResponsiveUnitDirection = ctc.ResponsiveUnitDirection;
    temp_ctc.ResponsiveStrength      = ctc.ResponsiveStrength;
    temp_ctc.CounterexampleUnitIdx   = ctc.CounterexampleUnitIdx;

    try
        temp_ctc.generateTrainLabels();
        temp_ctc.classifyUnits();

        % Assemble labels and embeddings from supervised UMAP
        all_reduction    = [temp_ctc.Reduction.Train; temp_ctc.Reduction.Test];
        train_labels_vec = temp_ctc.TrainLabels.sorted_y_train';
        test_labels_vec  = temp_ctc.UnitLabels(temp_ctc.TrainLabels.umap_test_idx)';
        all_labels       = [train_labels_vec; test_labels_vec];

        valid = ~isnan(all_labels);
        if sum(valid) < 10 || numel(unique(all_labels(valid))) < 2
            loss = 1;
            return
        end

        % Compute selected metric
        switch p_bo.ObjectiveMetric
            case "qf_dissimilarity"
                extras = temp_ctc.Reduction.Extras;
                if isempty(extras) || isempty(extras.qfd)
                    loss = 1;
                    return
                end
                [~, avgOverlap] = extras.getMatchSummary(4);
                if isnan(avgOverlap)
                    loss = 1;
                    return
                end
                loss = 1 - avgOverlap / 100;

            case "silhouette"
                sil_vals = silhouette(all_reduction(valid, :), all_labels(valid));
                loss = -mean(sil_vals);

            case "combined"
                sil_vals = silhouette(all_reduction(valid, :), all_labels(valid));
                extras = temp_ctc.Reduction.Extras;
                if ~isempty(extras) && ~isempty(extras.qfd)
                    [~, avgOverlap] = extras.getMatchSummary(4);
                else
                    avgOverlap = 0;
                end
                if isnan(avgOverlap); avgOverlap = 0; end
                loss = -mean(sil_vals) + (1 - avgOverlap / 100);

            otherwise
                error('CellTypeClassifier:unknownMetric', ...
                    'Unknown ObjectiveMetric: "%s".', p_bo.ObjectiveMetric);
        end

        % Optional interneuron fraction penalty
        inh_frac = sum(all_labels(valid) == 2) / sum(valid);
        if p_bo.UseInterneuronPenalty
            loss = loss + p_bo.InterneuronWeight * (inh_frac - p_bo.InterneuronTarget)^2;
        end

        if opts.Verbose
            active_vals = '';
            for vi = 1:numel(var_names)
                active_vals = [active_vals, sprintf('%s=', var_names{vi}), ...
                    sprintf('%.3g ', x.(var_names{vi}))]; %#ok<AGROW>
            end
            fprintf('  %s-> loss=%.3f inh=%.1f%%\n', active_vals, loss, 100*inh_frac);
        end
    catch ME
        warning('CellTypeClassifier:optimizeHyperparameters', ...
            'Evaluation failed: %s', ME.message);
        loss = 1;
    end
end

% -- Run Bayesian optimization ------------------------------------------------
if opts.Verbose
    fprintf('Bayesian optimization: %d topology variables (%s), %d evaluations, metric=%s\n', ...
        n_vars, strjoin(var_names, ', '), max_evals, p_bo.ObjectiveMetric);
end

bo_result = bayesopt(@objective, vars, ...
    'MaxObjectiveEvaluations', max_evals, ...
    'IsObjectiveDeterministic', false, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'Verbose', double(opts.Verbose), ...
    'PlotFcn', {});

% -- Extract best parameters -------------------------------------------------
best = bo_result.XAtMinObjective;

best_umap = struct('MinDist', best.MinDist, 'Spread', best.Spread);
if ismember('NNeighbors',          var_names); best_umap.NNeighbors          = best.NNeighbors; end
if ismember('SupervisedNNeighbors',var_names); best_umap.SupervisedNNeighbors = best.SupervisedNNeighbors; end
if ismember('NDims',               var_names); best_umap.NDims               = best.NDims; end

results.bestParams     = struct('UMAP', best_umap);
results.bestObjective  = bo_result.MinObjective;
results.bayesoptResult = bo_result;
results.allObjectives  = bo_result.XTrace;
results.nVars          = n_vars;

if opts.Verbose
    fprintf('\nBest topology parameters (metric=%s, %d vars optimized):\n', ...
        p_bo.ObjectiveMetric, n_vars);
    for vi = 1:numel(var_names)
        fprintf('  %-25s %g\n', var_names{vi}, best.(var_names{vi}));
    end
    fprintf('  Best objective:         %.4f\n', results.bestObjective);
end
end
