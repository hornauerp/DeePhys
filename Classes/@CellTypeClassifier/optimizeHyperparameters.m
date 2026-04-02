function results = optimizeHyperparameters(ctc, opts)
% OPTIMIZEHYPERPARAMETERS  Bayesian optimization of UMAP and outlier detection parameters.
%
% Runs bayesopt to find UMAP and outlier detection hyperparameters that
% maximize cluster quality using a selectable objective metric.
%
% Requires identifyResponsiveUnits() to have been run first (ground truth labels
% are fixed during optimization — only the UMAP/clustering parameters are tuned).
%
% When Auto* flags are enabled, those parameters are excluded from the Bayesian
% search, focusing optimization on the remaining free parameters. MaxEvaluations
% is adjusted to max(15, 10 * n_remaining_vars).
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
% OBJECTIVE METRICS (set via Parameters.BayesianOptimization.ObjectiveMetric):
%   "silhouette"       — negative mean silhouette on supervised 2D embedding
%   "qf_dissimilarity" — 1 - QF overlap from UMAP toolbox
%   "combined"         — negative silhouette + (1 - QF overlap/100)
%
% OPTIMIZED VARIABLES (up to 8, conditionally excluded by Auto* flags):
%   NDims                  — unsupervised UMAP output dimensions (excluded if AutoNDims)
%   NNeighbors             — unsupervised UMAP n_neighbors (excluded if AutoNNeighbors)
%   MinDist                — UMAP min_dist (always included)
%   Spread                 — UMAP spread (always included)
%   SupervisedNNeighbors   — supervised UMAP n_neighbors (excluded if AutoNNeighbors)
%   TargetWeight           — supervised UMAP target weight (excluded if AutoTargetWeight)
%   ContaminationFraction  — outlier detection contamination (always included)
%   CounterexampleRatio    — training ratio (excluded if AutoCounterexampleRatio)
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
p_outlr = ctc.Parameters.OutlierDetection;

% -- Define optimizable variables (conditionally, based on Auto* flags) ------
vars = optimizableVariable.empty;
var_names = {};

if ~p_umap.AutoNDims
    vars(end+1) = optimizableVariable('NDims', p_bo.NDimsRange, 'Type', 'integer');
    var_names{end+1} = 'NDims';
end
if ~p_umap.AutoNNeighbors
    vars(end+1) = optimizableVariable('NNeighbors', p_bo.NNeighborsRange, 'Type', 'integer');
    var_names{end+1} = 'NNeighbors';
end
% MinDist and Spread always included (no auto version)
vars(end+1) = optimizableVariable('MinDist', p_bo.MinDistRange, 'Transform', 'log');
var_names{end+1} = 'MinDist';
vars(end+1) = optimizableVariable('Spread', p_bo.SpreadRange);
var_names{end+1} = 'Spread';

if ~p_umap.AutoNNeighbors
    vars(end+1) = optimizableVariable('SupervisedNNeighbors', p_bo.SupervisedNNeighborsRange, 'Type', 'integer');
    var_names{end+1} = 'SupervisedNNeighbors';
end
if ~p_umap.AutoTargetWeight
    vars(end+1) = optimizableVariable('TargetWeight', p_bo.TargetWeightRange);
    var_names{end+1} = 'TargetWeight';
end
% ContaminationFraction always included (no auto version)
vars(end+1) = optimizableVariable('ContaminationFraction', p_bo.ContaminationRange);
var_names{end+1} = 'ContaminationFraction';

if ~p_outlr.AutoCounterexampleRatio
    vars(end+1) = optimizableVariable('CounterexampleRatio', p_bo.CounterexampleRatioRange, 'Type', 'integer');
    var_names{end+1} = 'CounterexampleRatio';
end

n_vars    = numel(vars);
max_evals = max(15, 10 * n_vars);

% -- Objective function -------------------------------------------------------
function loss = objective(x)
    temp_params = ctc.Parameters;

    % Apply candidate values for each active variable
    if ismember('NDims',               var_names); temp_params.UMAP.NDims               = x.NDims; end
    if ismember('NNeighbors',          var_names); temp_params.UMAP.NNeighbors          = x.NNeighbors; end
    if ismember('MinDist',             var_names); temp_params.UMAP.MinDist             = x.MinDist; end
    if ismember('Spread',              var_names); temp_params.UMAP.Spread              = x.Spread; end
    if ismember('SupervisedNNeighbors',var_names); temp_params.UMAP.SupervisedNNeighbors = x.SupervisedNNeighbors; end
    if ismember('TargetWeight',        var_names); temp_params.UMAP.TargetWeight        = x.TargetWeight; end
    if ismember('ContaminationFraction',var_names); temp_params.OutlierDetection.ContaminationFraction = x.ContaminationFraction; end
    if ismember('CounterexampleRatio', var_names); temp_params.OutlierDetection.CounterexampleRatio   = x.CounterexampleRatio; end

    temp_ctc = CellTypeClassifier(ctc.FeatureStore, ctc.UnitDataArray, temp_params);
    temp_ctc.ResponsiveUnitIdx = ctc.ResponsiveUnitIdx;

    try
        temp_ctc.generateTrainLabels();
        temp_ctc.classifyUnits();

        % Assemble labels and embeddings
        all_reduction = [temp_ctc.Reduction.Train; temp_ctc.Reduction.Test];
        train_labels_vec = temp_ctc.TrainLabels.sorted_y_train';
        test_labels_vec  = temp_ctc.UnitLabels(temp_ctc.TrainLabels.umap_test_idx)';
        all_labels = [train_labels_vec; test_labels_vec];

        valid = ~isnan(all_labels);
        if sum(valid) < 10 || numel(unique(all_labels(valid))) < 2
            loss = 1;
            return
        end

        % Compute selected metric
        switch p_bo.ObjectiveMetric
            case "silhouette"
                sil_vals = silhouette(all_reduction(valid, :), all_labels(valid));
                loss = -mean(sil_vals);

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
    fprintf('Bayesian optimization: %d variables (%s), %d evaluations, metric=%s\n', ...
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

best_umap   = struct('MinDist', best.MinDist, 'Spread', best.Spread);
best_outlr  = struct('ContaminationFraction', best.ContaminationFraction);

if ismember('NDims',               var_names); best_umap.NDims               = best.NDims; end
if ismember('NNeighbors',          var_names); best_umap.NNeighbors          = best.NNeighbors; end
if ismember('SupervisedNNeighbors',var_names); best_umap.SupervisedNNeighbors = best.SupervisedNNeighbors; end
if ismember('TargetWeight',        var_names); best_umap.TargetWeight         = best.TargetWeight; end
if ismember('CounterexampleRatio', var_names); best_outlr.CounterexampleRatio = best.CounterexampleRatio; end

results.bestParams     = struct('UMAP', best_umap, 'OutlierDetection', best_outlr);
results.bestObjective  = bo_result.MinObjective;
results.bayesoptResult = bo_result;
results.allObjectives  = bo_result.XTrace;
results.nVars          = n_vars;

if opts.Verbose
    fprintf('\nBest parameters (metric=%s, %d vars optimized):\n', ...
        p_bo.ObjectiveMetric, n_vars);
    for vi = 1:numel(var_names)
        fprintf('  %-25s %g\n', var_names{vi}, best.(var_names{vi}));
    end
    fprintf('  Best objective:         %.4f\n', results.bestObjective);
end
end
