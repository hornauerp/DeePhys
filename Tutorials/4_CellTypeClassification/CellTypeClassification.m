%% Tutorial 4 — Cell Type Classification
%
% Classifies neurons as excitatory (1) or inhibitory (2) using a supervised
% UMAP pipeline built around two ground-truth strategies:
%
%   Drug-response (default): Bootstrap firing-rate test identifies units that
%   increase firing after stimulus — these are the inhibitory ground-truth set.
%   Counterexamples are selected farthest from the inhibitory centroid in feature space.
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
% GroundTruthMethod controls how inhibitory candidates are identified:
%   'two_window' (default): bootstrap permutation test comparing pre vs post
%     firing rate within each culture. Requires distinct baseline and treatment
%     recordings selectable by GroupingVar.
%   'full_curve': per-unit Spearman rank correlation between firing rate and
%     GroupingVar across all recordings of the culture. Requires the same
%     UnitIDs across recordings. More robust for multi-concentration dose-response
%     experiments (e.g. 0, 1, 10, 100 µM). Falls back to 'two_window' if fewer
%     than Bootstrap.MinRecordings recordings are available.
%   'metadata': read labels directly from UnitTable — see Section 6.
params.Bootstrap.GroundTruthMethod = 'two_window';

% Alpha: 1e-10 is deliberately strict. The downstream PCA outlier detection
%   handles impure candidates, but over-inclusion here (weak alpha) floods the
%   training set with ambiguous units and degrades embedding quality more than
%   under-inclusion. If you get fewer than ~10 responsive units per culture,
%   consider relaxing to 1e-6 before enabling UseFDR.
%
% FullCurveAlpha: equivalent threshold for 'full_curve' (Spearman p-value).
%   0.05 is appropriate since the Spearman test is already conservative for
%   monotonic dose-response and the small N (usually 4-6 recordings).
%
% UseFDR: Not recommended as default. BH correction scales with test count and
%   can be over-conservative on small datasets. The fixed alpha with strict
%   threshold is simpler and better calibrated to this specific use case.
%
% Direction: 'increase' targets inhibitory units (disinhibition → FR increase).
%   Use 'both' only if your stimulus elicits both excitatory and inhibitory
%   activity-dependent responses simultaneously (uncommon in acute pharmacology).
params.Bootstrap.Alpha          = 1e-10;
params.Bootstrap.FullCurveAlpha = 0.05;
params.Bootstrap.NIter          = 1000;
params.Bootstrap.Direction      = 'increase';

% ── Normalization and recording selection ──────────────────────────────────────
%
% Per-chip z-score before the global z-score removes chip-level offsets in
% feature distributions (different electrode impedances, culture densities).
% Set NormalizationVar = '' to skip per-group normalization on single-chip datasets.
params.UMAP.NormalizationVar = 'ChipID';

% GroupingVar / GroupingValues: controls which recordings contribute to the
% feature matrix and the bootstrap test. Only recordings whose GroupingVar
% value is in GroupingValues are included. For dose-response experiments,
% GroupingValues = 0 (baseline/vehicle) is correct — you want features from
% the untreated state so that drug-response changes do not contaminate the
% feature space. If you want features from all recordings, set GroupingValues
% to [] or to the full range. Must match the metadata column name exactly.
params.UMAP.GroupingVar    = 'Concentration';
params.UMAP.GroupingValues = 0;

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
params.UMAP.MinDist = 0.1;
params.UMAP.Spread  = 1.0;

% ── UMAP neighbor counts ───────────────────────────────────────────────────────
%
% AutoNNeighbors = true: sets n_neighbors for both the unsupervised and supervised
%   UMAP passes to max(MinNNeighbors, sqrt(N)), scaling automatically with dataset
%   size. Recommended — prevents n_neighbors from exceeding N on small datasets or
%   being too local on large ones. When enabled, explicit NNeighbors is ignored.
%
% SupervisedNDims = 2: embedding dimensionality for supervised UMAP. 2 gives a
%   human-interpretable scatter plot. BayOpt can explore higher values (up to 10).
%
% SupervisedNNeighbors is additionally clamped at runtime by a minority-class
%   ceiling: floor(n_inh / 5). This prevents manifold tearing when the inhibitory
%   training set is small. BayOpt searches the range [5, 50]; the ceiling is then
%   applied to whatever value it selects.
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
% 0.3 (default) is deliberately conservative — drug-response labels are noisy
% (network effects, spike sorting errors, contaminated responsive units) and a
% lower weight prevents the embedding from forcing false separation.
% BayOpt can tune this jointly with the topology parameters (Section 8).
params.UMAP.TargetWeight = 0.3;

% ── Counterexample selection ───────────────────────────────────────────────────
%
% CounterexampleRatio = 1: one excitatory candidate per inhibitory candidate.
% Counterexamples are selected by farthest-from-inhibitory-centroid in normalized
% feature space (correlation distance). This ensures the excitatory training set
% is well-separated from the inhibitory cluster and avoids boundary-ambiguous units.
%
% Increase ratio to 2–3 if you consistently observe too many false inhibitory
% labels (inhibitory fraction > 30%).
params.OutlierDetection.CounterexampleRatio = 1;

% ── Outlier detection domain ──────────────────────────────────────────────────
%
% Outlier detection runs in PCA space (top 15 components, default) rather than
% UMAP space. This decouples training label generation from UMAP hyperparameters,
% avoiding a circular dependency where BayOpt changes indirectly affect labels.
% Set Domain = 'umap' to use the unsupervised UMAP embedding instead.
params.OutlierDetection.Domain = 'pca';

% ── Reproducibility ───────────────────────────────────────────────────────────
%
% RNGSeed controls the random seed used for UMAP and bootstrap permutations.
% Change to verify stability: if results change substantially across seeds,
% either training set is too small or the embedding geometry is degenerate
% (in which case run optimizeHyperparameters). Default: 42.
params.RNGSeed = 42;

% ── Confidence threshold (optional) ───────────────────────────────────────────
%
% By default all units are classified even if the kNN vote is narrow.
% Enable UseConfidenceThreshold to mark low-confidence predictions as NaN.
% ConfidenceThreshold = 0.3 is a conservative floor (near-random boundary for
% binary classification). Units below this threshold are genuinely ambiguous.
% Adjust after inspecting ctc.UnitConfidence distribution (Section 9).
%
% params.Classification.UseConfidenceThreshold = true;
% params.Classification.ConfidenceThreshold    = 0.3;

% ── Ensemble classification (optional) ────────────────────────────────────────
%
% Runs the full pipeline N times with different UMAP seeds and assigns labels
% by majority vote. Confidence = fraction of seeds agreeing on the winning class.
% MinAgreement sets the minimum vote fraction for label assignment — units below
% it are marked NaN. Substantially increases runtime (Nx pipeline cost).
%
% Use when single-run results are unstable (ARI < 0.9 in assessStability) or
% when downstream analysis is sensitive to occasional label flips.
%
% params.Ensemble.Enabled      = true;
% params.Ensemble.Seeds        = [42, 1042, 2042, 3042, 4042];
% params.Ensemble.MinAgreement = 0.6;

% ── Diagnostics (optional) ────────────────────────────────────────────────────
%
% When enabled, each pipeline stage automatically generates a diagnostic figure:
%   identifyResponsiveUnits → diagnosticResponsiveUnits  (responsive unit check)
%   generateTrainLabels     → diagnosticTrainLabels       (label quality check)
%   classifyUnits           → diagnosticClassification    (prediction plausibility)
%   optimizeHyperparameters → diagnosticOptimization      (BayOpt convergence)
%
% Each diagnostic can also be called manually: ctc.diagnosticResponsiveUnits()
%
% params.Diagnostics.Enable      = true;    % master toggle (default false)
% params.Diagnostics.SaveDir     = '/path/to/output';   % save PNGs here
% params.Diagnostics.ShowFigures = true;    % false suppresses display (batch mode)

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
% Internally calls buildNormalizedFeatures (shared with classifyUnits):
%   - Extracts harmonized waveforms + ACGs once; builds feature matrix
%   - Applies per-group z-score (NormalizationVar groups)
%   - HarmonizedWaveforms/HarmonizedACGs are available after this call
%
% Then on the training subset:
%   1. Fits global z-score + scaling; stores NormalizationParams
%   2. Optional feature selection (FeatureSelection flag)
%   3. Unsupervised UMAP embedding (used for outlier detection)
%   4. Detects and removes outliers from responsive candidates (PCA space by default)
%   5. Selects counterexamples farthest from inhibitory centroid in feature space
%   6. Detects and removes outliers from counterexamples (PCA space by default)
%      with top-up from a reserve pool to maintain the target count

ctc.generateTrainLabels();

tl = ctc.TrainLabels;
fprintf('Train set: %d excitatory, %d inhibitory\n', ...
    sum(tl.sorted_y_train == 1), sum(tl.sorted_y_train == 2));

%% 5  Classify all units (supervised UMAP)
%
% Reuses the NormalizedFeatures cache from generateTrainLabels (no second
% feature extraction). Applies the stored NormalizationParams (global
% z-score, NaN-column removal, scaling) to all unique units, then trains
% the supervised UMAP and projects test units into the embedding.
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
%
% HarmonizedWaveforms and HarmonizedACGs are populated by generateTrainLabels
% (via the shared feature cache). They are already available after Section 4.

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
% overlap completely. These artifacts arise from UMAP hyperparameters and feature
% group imbalances interacting badly with a particular dataset.
%
% What is optimized: 7 parameters (MinDist, Spread, SupervisedNDims,
%   SupervisedNNeighbors, TargetWeight, ACGWeight, WaveformWeight).
%
%   ACGWeight / WaveformWeight: the ACG feature group has ~1001 dimensions vs
%   ~241 for waveforms. Without correction UMAP treats each dimension equally,
%   so ACG dominates ~81% of all pairwise distances. BayOpt finds the optimal
%   per-group scaling; starting values default to 1.0 (equal group contribution).
%   These cannot be chosen by intuition — only BayOpt can find a non-trivial value.
%
%   SupervisedNNeighbors is additionally clamped at runtime by the minority
%   class ceiling (floor(n_inh / 5)) to prevent manifold tearing.
%
% Objective metric: weighted silhouette (default) — evaluates cluster separation
%   on the test embedding only, with asymmetric weighting (0.7 inhibitory +
%   0.3 excitatory). The inhibitory class is harder to classify (minority class),
%   so it receives more weight. Evaluated on test units only to avoid conflating
%   parameter quality with the strength of the supervised signal in train units.
%   Trustworthiness, silhouette, and combined are also available.
%
% Run count: 60 evaluations (10 per variable).
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
