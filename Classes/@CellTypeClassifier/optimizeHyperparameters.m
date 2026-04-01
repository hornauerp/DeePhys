function results = optimizeHyperparameters(ctc, opts)
% OPTIMIZEHYPERPARAMETERS  Bayesian optimization of UMAP and outlier detection parameters.
%
% Runs bayesopt to find UMAP and outlier detection hyperparameters that
% maximize cluster quality using a selectable objective metric.
%
% Requires identifyResponsiveUnits() to have been run first (ground truth labels
% are fixed during optimization — only the UMAP/clustering parameters are tuned).
%
% USAGE:
%   ctc = CellTypeClassifier(rg, params);
%   ctc.identifyResponsiveUnits();
%   results = ctc.optimizeHyperparameters();
%   % Apply best parameters:
%   ctc.Parameters = parseStructParameters(ctc.Parameters, results.bestParams);
%   ctc.generateTrainLabels();
%   ctc.classifyUnits();
%
% NAME-VALUE:
%   Verbose - Print progress during optimization (default true)
%
% OBJECTIVE METRICS (set via Parameters.BayesianOptimization.ObjectiveMetric):
%   "silhouette"       — negative mean silhouette on supervised 2D embedding (default)
%   "qf_dissimilarity" — 1 - QF overlap from UMAP toolbox (train-test cluster match)
%   "combined"         — negative silhouette + (1 - QF overlap/100)
%
% OPTIMIZED VARIABLES (8 total):
%   NDims                  — unsupervised UMAP output dimensions (affects isolation forest)
%   NNeighbors             — unsupervised UMAP n_neighbors
%   MinDist                — UMAP min_dist (shared unsupervised/supervised)
%   Spread                 — UMAP spread (affects cluster separation)
%   SupervisedNNeighbors   — supervised UMAP n_neighbors
%   TargetWeight           — supervised UMAP balance between data and label topology
%   ContaminationFraction  — isolation forest contamination
%   CounterexampleRatio    — excitatory:inhibitory training set ratio
%
% OUTPUTS:
%   results - struct with fields:
%     .bestParams     - struct of optimal parameter values (nested for direct merging)
%     .bestObjective  - best (lowest) objective value achieved
%     .bayesoptResult - full BayesianOptimization object for diagnostics
%     .allObjectives  - table of all evaluations

arguments
    ctc CellTypeClassifier
    opts.Verbose (1,1) logical = true
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before optimizeHyperparameters()');

p_bo = ctc.Parameters.BayesianOptimization;

% ── Define optimizable variables ─────────────────────────────────────────────
vars = [
    optimizableVariable('NDims',                 p_bo.NDimsRange,                 'Type', 'integer')
    optimizableVariable('NNeighbors',            p_bo.NNeighborsRange,            'Type', 'integer')
    optimizableVariable('MinDist',               p_bo.MinDistRange,               'Transform', 'log')
    optimizableVariable('Spread',                p_bo.SpreadRange)
    optimizableVariable('SupervisedNNeighbors',  p_bo.SupervisedNNeighborsRange,  'Type', 'integer')
    optimizableVariable('TargetWeight',          p_bo.TargetWeightRange)
    optimizableVariable('ContaminationFraction', p_bo.ContaminationRange)
    optimizableVariable('CounterexampleRatio',   p_bo.CounterexampleRatioRange,   'Type', 'integer')
];

% ── Objective function ───────────────────────────────────────────────────────
function loss = objective(x)
    % Create a temporary CTC with candidate parameters
    temp_params = ctc.Parameters;
    temp_params.UMAP.NDims                 = x.NDims;
    temp_params.UMAP.NNeighbors            = x.NNeighbors;
    temp_params.UMAP.MinDist               = x.MinDist;
    temp_params.UMAP.Spread                = x.Spread;
    temp_params.UMAP.SupervisedNNeighbors  = x.SupervisedNNeighbors;
    temp_params.UMAP.TargetWeight          = x.TargetWeight;
    temp_params.OutlierDetection.ContaminationFraction = x.ContaminationFraction;
    temp_params.OutlierDetection.CounterexampleRatio   = x.CounterexampleRatio;

    temp_ctc = CellTypeClassifier(ctc.RecordingGroup, temp_params);
    temp_ctc.UnitList          = ctc.UnitList;
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
                    'Unknown ObjectiveMetric: "%s". Use "silhouette", "qf_dissimilarity", or "combined".', ...
                    p_bo.ObjectiveMetric);
        end

        % Optional interneuron fraction penalty
        inh_frac = sum(all_labels(valid) == 2) / sum(valid);
        if p_bo.UseInterneuronPenalty
            loss = loss + p_bo.InterneuronWeight * (inh_frac - p_bo.InterneuronTarget)^2;
        end

        if opts.Verbose
            fprintf('  NDims=%d NNeigh=%d MinDist=%.3f Spread=%.2f SupNNeigh=%d TgtWt=%.2f Contam=%.2f CERatio=%d → loss=%.3f inh=%.1f%%\n', ...
                x.NDims, x.NNeighbors, x.MinDist, x.Spread, ...
                x.SupervisedNNeighbors, x.TargetWeight, ...
                x.ContaminationFraction, x.CounterexampleRatio, ...
                loss, 100*inh_frac);
        end
    catch ME
        warning('CellTypeClassifier:optimizeHyperparameters', ...
            'Evaluation failed: %s', ME.message);
        loss = 1;
    end
end

% ── Run Bayesian optimization ────────────────────────────────────────────────
if opts.Verbose
    fprintf('Starting Bayesian optimization (%d evaluations, metric=%s, 8 variables)...\n', ...
        p_bo.MaxEvaluations, p_bo.ObjectiveMetric);
end

bo_result = bayesopt(@objective, vars, ...
    'MaxObjectiveEvaluations', p_bo.MaxEvaluations, ...
    'IsObjectiveDeterministic', false, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'Verbose', double(opts.Verbose), ...
    'PlotFcn', {});

% ── Extract best parameters (nested struct for parseStructParameters) ────────
best = bo_result.XAtMinObjective;
results.bestParams = struct( ...
    'UMAP', struct( ...
        'NDims',                 best.NDims, ...
        'NNeighbors',            best.NNeighbors, ...
        'MinDist',               best.MinDist, ...
        'Spread',                best.Spread, ...
        'SupervisedNNeighbors',  best.SupervisedNNeighbors, ...
        'TargetWeight',          best.TargetWeight), ...
    'OutlierDetection', struct( ...
        'ContaminationFraction', best.ContaminationFraction, ...
        'CounterexampleRatio',   best.CounterexampleRatio));
results.bestObjective  = bo_result.MinObjective;
results.bayesoptResult = bo_result;
results.allObjectives  = bo_result.XTrace;

if opts.Verbose
    fprintf('\nBest parameters found (metric=%s):\n', p_bo.ObjectiveMetric);
    fprintf('  NDims:                 %d\n', best.NDims);
    fprintf('  NNeighbors:            %d\n', best.NNeighbors);
    fprintf('  MinDist:               %.4f\n', best.MinDist);
    fprintf('  Spread:                %.3f\n', best.Spread);
    fprintf('  SupervisedNNeighbors:  %d\n', best.SupervisedNNeighbors);
    fprintf('  TargetWeight:          %.3f\n', best.TargetWeight);
    fprintf('  ContaminationFraction: %.3f\n', best.ContaminationFraction);
    fprintf('  CounterexampleRatio:   %d\n', best.CounterexampleRatio);
    fprintf('  Best objective:        %.4f\n', results.bestObjective);
end
end
