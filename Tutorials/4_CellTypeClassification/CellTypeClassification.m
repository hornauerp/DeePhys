%% Tutorial 4 — Cell Type Classification
%
% Classifies neurons as excitatory (1) or inhibitory (2) using a supervised
% UMAP pipeline built around two ground-truth strategies:
%
%   Drug-response (default): Bootstrap firing-rate test identifies units that
%   increase firing after stimulus — these are the inhibitory ground-truth set.
%   Counterexamples are sampled uniformly from the non-responsive pool.
%
%   Metadata (pure cultures): Both classes supplied directly from a UnitTable
%   column (e.g. optogenetics tag, genetic marker). Skips bootstrap entirely.
%
% Pipeline: identifyResponsiveUnits → generateTrainLabels → classifyUnits
%
% Prerequisites:
%   - Tutorial 1 completed (FeatureStore and RecordingProcessors saved)
%   - DeePhys on the MATLAB path

%% 1  Load data

fs_file    = '/path/to/FeatureStore.mat';
proc_paths = {'/path/to/proc1.mat', '/path/to/proc2.mat'};

fs    = FeatureStore.load(fs_file);
procs = RecordingProcessor.loadMany(proc_paths);

% Build UnitData array (must match UnitTable row order)
ud = [];
for i = 1:numel(procs)
    ud = [ud, procs(i).Units]; %#ok<AGROW>
end

%% 2  Recommended parameter set
%
% Parameters below reflect considered defaults for the typical MEA drug-response
% experiment. Each choice is annotated with what it controls and what to try if
% the result is unsatisfactory.

params = struct();

% ── Harmonization ─────────────────────────────────────────────────────────────
%
% ACGBinSize: 0.5 ms captures fine rebound inhibition peaks. Increase to 1 ms
%   for speed on large datasets. Do not go below 0.2 ms — bins become too sparse
%   for low-firing units.
% ACGLag: ±100 ms covers the full rebound window for most inhibitory subtypes.
%   Reduce to ±50 ms if memory is limited, but you lose late-rebound features.
% ACGSource: 'FullACG' uses the parent recording ACG (better statistics than
%   the per-recording ACG). Falls back to per-recording if Parent_ACG* is absent.
params.Harmonization.ACGBinSize = 0.0005;
params.Harmonization.ACGLag     = 0.1;
params.Harmonization.ACGSource  = 'FullACG';

% ── Bootstrap firing-rate test ─────────────────────────────────────────────────
%
% Alpha: 1e-10 is deliberately strict. The downstream UMAP outlier detection
%   handles impure candidates, but over-inclusion here (weak alpha) floods the
%   training set with ambiguous units and degrades embedding quality more than
%   under-inclusion. If you get fewer than ~10 responsive units per culture,
%   consider relaxing to 1e-6 before enabling UseFDR.
%
% UseFDR: Not recommended as default. BH correction scales with test count and
%   can be over-conservative on small datasets. The fixed alpha with strict
%   threshold is simpler and better calibrated to this specific use case.
%
% Direction: 'increase' targets inhibitory units (disinhibition → FR increase).
%   Use 'both' only if your stimulus elicits both excitatory and inhibitory
%   activity-dependent responses simultaneously (uncommon in acute pharmacology).
params.Bootstrap.Alpha     = 1e-10;
params.Bootstrap.NIter     = 1000;
params.Bootstrap.Direction = 'increase';

% ── Normalization ──────────────────────────────────────────────────────────────
%
% Per-chip z-score before the global z-score removes chip-level offsets in
% feature distributions (different electrode impedances, culture densities).
% Set NormalizationVar = '' to skip per-group normalization on single-chip datasets.
params.UMAP.NormalizationVar = 'ChipID';

% ── UMAP geometry ─────────────────────────────────────────────────────────────
%
% These parameters control embedding shape. The defaults below work well for
% typical datasets (100–2000 units). If the supervised embedding produces
% degenerate geometry (thin ellipses, collapsed clusters), run
% optimizeHyperparameters() — see Section 8.
%
% MinDist: 0.1 is a balanced default. Lower (0.05) packs clusters tighter and
%   can improve separation for clean datasets, but amplifies noise artifacts.
%   Higher (0.3) produces more spread-out embeddings that look better but may
%   blur the E/I boundary.
%
% Spread: 1.0 is the UMAP default. Adjust jointly with MinDist — changing one
%   without the other often produces no useful change.
%
% NNeighbors: 30 for the unsupervised embedding. More neighbors → more global
%   structure preserved; fewer → more local clusters. AutoNNeighbors = true
%   sets this to max(15, sqrt(N)) which is useful when dataset size varies
%   greatly across experiments.
params.UMAP.MinDist    = 0.1;
params.UMAP.Spread     = 1.0;
params.UMAP.NNeighbors = 30;

% ── Supervised UMAP dimensionality ────────────────────────────────────────────
%
% AutoSupervisedNDims uses the TWO-NN intrinsic dimensionality estimator
% (Facco et al. 2017) to set the number of supervised embedding dimensions.
% This is strongly preferred over fixing SupervisedNDims = 2:
%   - Real E/I separation often needs 3–6 dimensions to capture ACG, waveform,
%     and firing-rate features simultaneously.
%   - Forcing 2D compresses information the classifier needs and increases
%     sensitivity to seed and geometry artifacts.
%   - TWO-NN estimates d from the ratio of 2nd to 1st nearest-neighbor distances
%     under a Pareto null — no tuning required.
%
% Alternative: AutoSupervisedNDims = false, SupervisedNDims = 2.
%   Use this only when you need a human-interpretable 2D scatter plot and
%   are willing to accept lower classification accuracy.
%
% MaxSupervisedNDims: caps the estimate at 10 to prevent runaway on noisy data.
%   Rarely activates in practice (typical estimate is 3–6).
params.UMAP.AutoSupervisedNDims = true;
params.UMAP.MinSupervisedNDims  = 2;
params.UMAP.MaxSupervisedNDims  = 10;

% ── Supervised NNeighbors ─────────────────────────────────────────────────────
%
% AutoNNeighbors = true also applies to the supervised pass, setting
% SupervisedNNeighbors = max(15, sqrt(N_train)). This prevents k > N_train
% on small training sets.
%
% We recommend keeping this on: supervised UMAP training set size varies
% substantially across experiments and the sqrt heuristic is well-calibrated.
params.UMAP.AutoNNeighbors   = true;
params.UMAP.MinNNeighbors    = 15;

% ── kNN confidence ─────────────────────────────────────────────────────────────
%
% Distance-weighted kNN in UMAP space. AutoConfidenceK = true sets k = max(5,
% sqrt(N_train)), preventing k from approaching N_train (which collapses all
% confidences to the class prior).
%
% Alternative: AutoConfidenceK = false, ConfidenceK = 15.
%   Fine for large datasets; breaks on small ones.
params.UMAP.AutoConfidenceK = true;

% ── TargetWeight ──────────────────────────────────────────────────────────────
%
% Controls how strongly label information pulls the supervised topology.
% 0.5 (default) balances label influence against data geometry — appropriate
% for noisy drug-response labels where the training set is not pure.
%
% MetadataTargetWeight (0.8) is used automatically when ground truth comes from
% metadata (Section 6) because those labels are externally verified and cleaner.
%
% Do not set TargetWeight very high (>0.8) for drug-response labels: if labels
% are noisy the embedding will force false separation.
params.UMAP.TargetWeight         = 0.5;
params.UMAP.MetadataTargetWeight = 0.8;

% ── Counterexample selection ───────────────────────────────────────────────────
%
% CounterexampleRatio = 1: one excitatory candidate per inhibitory candidate.
% Uniform random selection from the non-responsive pool (not farthest-point,
% not boundary-weighted). Outlier detection in UMAP space cleans both classes
% after selection, so the sampling strategy only needs to be unbiased.
%
% Farthest-point sampling was rejected: it preferentially picks extreme
% feature values and artifact-like units that happen to be far from the
% responsive cluster but are not representative excitatory neurons.
%
% Increase ratio to 2–3 if you consistently observe too many false inhibitory
% labels (inhibitory fraction > 30%).
params.OutlierDetection.CounterexampleRatio = 1;

% Construct classifier with recommended parameters
ctc = CellTypeClassifier(fs, ud, params);

%% 3  Identify responsive units (drug-response scenario)
%
% Bootstrap permutation test: compares pre-stimulus vs post-stimulus firing
% rate across cultures. Units with a significant rate increase become the
% inhibitory candidate set (positive class).
%
% The filter argument restricts which cultures contribute candidates.
% Useful if only some conditions contain a drug that produces the response.

ctc.identifyResponsiveUnits();

% With filter — only use control wells (concentration = 0):
%   ctc.identifyResponsiveUnits({'Concentration', 0});

fprintf('Inhibitory candidates: %d / %d total units\n', ...
    sum(ctc.ResponsiveUnitIdx), numel(ctc.ResponsiveUnitIdx));

%% 4  Generate training labels
%
% Embeds all units with unsupervised UMAP, then:
%   1. Detects and removes outliers from responsive candidates (in UMAP space)
%   2. Selects counterexamples uniformly from the non-responsive pool
%   3. Detects and removes outliers from counterexamples (in UMAP space)
%      with top-up from a reserve pool to maintain the target count
%   4. Stores NormalizationParams for consistent application in classifyUnits
%
% SupervisedNDims is estimated here via TWO-NN if AutoSupervisedNDims = true.
% The estimate is logged to the console.

ctc.generateTrainLabels();

tl = ctc.TrainLabels;
fprintf('Train set: %d excitatory, %d inhibitory\n', ...
    sum(tl.sorted_y_train == 1), sum(tl.sorted_y_train == 2));

%% 5  Classify all units (supervised UMAP)
%
% Trains supervised UMAP on the labelled training set, projects all remaining
% units into the embedding, assigns labels by nearest-neighbour lookup.
% Confidence = distance-weighted kNN vote fraction in [0, 1].

ctc.classifyUnits();

labels = ctc.UnitLabels;   % 1=excitatory, 2=inhibitory, NaN=unclassified
n_exc  = sum(labels == 1, 'omitnan');
n_inh  = sum(labels == 2, 'omitnan');
n_nan  = sum(isnan(labels));
fprintf('Excitatory: %d  Inhibitory: %d  Unclassified: %d\n', n_exc, n_inh, n_nan);

% Typical inhibitory fraction in MEA cultures: 15–25%. Values outside this
% range suggest parameter issues — check the UMAP embedding in Section 9.

%% 6  Metadata-driven labels (pure culture scenario)
%
% When ground truth is externally supplied (optogenetics, genetic markers,
% morphology tags), use GroundTruthMethod = "metadata" to read labels
% directly from the UnitTable. This skips the bootstrap entirely.
%
% LabelField:           column in fs.UnitTable containing cell type labels
% ResponsiveClassValue: the string value that denotes the inhibitory class
%
% The MetadataTargetWeight (0.8) is applied automatically for this path
% because metadata labels are cleaner and warrant stronger label influence
% in the supervised UMAP.

params_meta = params;  % inherit all recommended parameters above
params_meta.Bootstrap.GroundTruthMethod    = 'metadata';
params_meta.Bootstrap.LabelField           = 'CellType';   % column in UnitTable
params_meta.Bootstrap.ResponsiveClassValue = 'inhibitory'; % string for inh class

ctc_meta = CellTypeClassifier(fs, ud, params_meta);
ctc_meta.identifyResponsiveUnits();   % reads labels from UnitTable, returns immediately
ctc_meta.generateTrainLabels();
ctc_meta.classifyUnits();

fprintf('Metadata path: %d exc, %d inh\n', ...
    sum(ctc_meta.UnitLabels==1,'omitnan'), ...
    sum(ctc_meta.UnitLabels==2,'omitnan'));

%% 7  Inspect harmonized features

wf  = ctc.HarmonizedWaveforms;   % (N_samples x N_units)
acg = ctc.HarmonizedACGs;         % (N_bins   x N_units)
sr  = ctc.HarmonizedSR;           % target waveform sampling rate

fprintf('Waveform : %d x %d at %.0f Hz\n', size(wf,1), size(wf,2), sr);
fprintf('ACG      : %d x %d\n', size(acg,1), size(acg,2));

plotCellTypeFeatures(ctc);
sortACGsByPeak(ctc);

%% 8  Bayesian optimization — fix degenerate UMAP geometry
%
% Use this when the supervised embedding produces poor geometry: thin ellipses
% along a single axis, collapsed clusters, or a layout where E/I populations
% overlap completely. These artifacts arise from MinDist/Spread/NNeighbors
% interacting badly with a particular dataset's feature density.
%
% What is optimized: topology parameters only (MinDist, Spread, NNeighbors,
%   SupervisedNNeighbors). Parameters that control label assignment
%   (TargetWeight, CounterexampleRatio) are excluded — they can trivially
%   improve the metric (tighter clusters) without improving classification.
%
% Objective metric: qf_dissimilarity (default) — measures topology-label
%   agreement via QF overlap from the UMAP toolbox. Not gameable by topology
%   parameters because it compares the embedding structure to the supplied
%   labels, not just cluster compactness. Silhouette is the alternative if
%   the QF metric is unavailable, but it is gameable by TargetWeight.
%
% Run count: max(15, 10 * n_vars) evaluations. With AutoNNeighbors = true,
%   NNeighbors and SupervisedNNeighbors are fixed, leaving 2 free variables
%   (MinDist, Spread) → 20 evaluations. With AutoNNeighbors = false, 4 free
%   variables → 40 evaluations.
%
% This should be treated as a one-time calibration step for a new experimental
% preparation, not run on every dataset.

ctc_opt = CellTypeClassifier(fs, ud, params);
ctc_opt.identifyResponsiveUnits();

results = ctc_opt.optimizeHyperparameters();
fprintf('Optimized %d variables, best objective: %.4f\n', ...
    results.nVars, results.bestObjective);

% Apply best topology parameters, then run the rest of the pipeline
ctc_opt.Parameters = parseStructParameters(ctc_opt.Parameters, results.bestParams);
ctc_opt.generateTrainLabels();
ctc_opt.classifyUnits();

% Check stability with the optimized parameters
stability = ctc_opt.assessStability('NRuns', 5);
fprintf('Stability: ARI = %.3f +/- %.3f\n', stability.meanARI, stability.stdARI);

%% 9  Validate training labels
%
% Leave-one-culture-out cross-validation on the training set. Trains on all
% cultures except one, classifies held-out units, and reports accuracy.
% Low per-culture accuracy (<75%) indicates that culture's training labels
% are inconsistent with the rest of the dataset — check identifyResponsiveUnits
% results for that culture.

ctc.validateTrainingLabels();

%% 10  Inspect UMAP embedding
%
% The supervised embedding is stored in ctc.Reduction.Train (training units)
% and ctc.Reduction.Test (all other units). Colour by label to verify separation.

figure;
subplot(1,2,1);
scatter(ctc.Reduction.Train(:,1), ctc.Reduction.Train(:,2), 10, ...
    ctc.TrainLabels.sorted_y_train, 'filled');
colormap([0.2 0.6 0.9; 0.9 0.3 0.2]);
title('Training embedding'); xlabel('UMAP 1'); ylabel('UMAP 2');
colorbar('Ticks', [1 2], 'TickLabels', {'Exc', 'Inh'});

subplot(1,2,2);
test_embed = ctc.Reduction.Test;
test_label = ctc.UnitLabels(ctc.TrainLabels.umap_test_idx);
scatter(test_embed(:,1), test_embed(:,2), 10, test_label, 'filled');
colormap([0.2 0.6 0.9; 0.9 0.3 0.2]);
title('Test embedding'); xlabel('UMAP 1');
colorbar('Ticks', [1 2], 'TickLabels', {'Exc', 'Inh'});

% Signs of good geometry:
%   - Two visually distinct clusters along the first UMAP axis
%   - Training units form compact cores within each cluster
%   - No thin ellipses or single-axis layout (if present, run Section 8)
%
% Confidence threshold: if many test units have low confidence (<0.5) and are
% scattered throughout the embedding, consider increasing Bootstrap.Alpha
% (stricter responsive threshold) to improve training set purity.
figure;
scatter(test_embed(:,1), test_embed(:,2), 10, ...
    ctc.UnitConfidence(ctc.TrainLabels.umap_test_idx), 'filled');
colormap(parula); colorbar;
title('Test confidence'); xlabel('UMAP 1'); ylabel('UMAP 2');
