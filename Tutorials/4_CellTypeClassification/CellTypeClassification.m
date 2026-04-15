%% Tutorial 4 — Cell Type Classification
%
% Classifies neurons as excitatory (1) or inhibitory (2) using a transductive
% graph-based pipeline. All units are embedded together in a single unsupervised
% UMAP; Louvain community detection identifies inhibitory communities from their
% responsive-unit enrichment; label propagation on the UMAP graph assigns cell
% types to every unit without a separate test-set projection.
%
%   Drug-response (default): Bootstrap firing-rate test identifies units that
%   increase firing after stimulus — these are the inhibitory ground-truth set.
%   Louvain communities are evaluated by their fraction of responsive units.
%
%   Metadata (pure cultures): Labels read directly from a UnitTable column
%   (e.g. optogenetics tag, genetic marker). Skips bootstrap entirely.
%
% Pipeline:
%   identifyResponsiveUnits()       ← flag drug-responsive (inhibitory) candidates
%   [optimizeUnsupervisedUMAP()]    ← optional Phase 1 BayOpt (recommended)
%   generateTrainLabels()           ← Louvain on UMAP graph → inh/CE training set
%   classifyUnits()                 ← label propagation on UMAP graph (transductive)
%
% Prerequisites:
%   - Tutorial 1 completed (FeatureStore and RecordingProcessors saved)
%   - DeePhys on the MATLAB path
%   - Brain Connectivity Toolbox (BCT) on the MATLAB path (community_louvain)

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
%
% A minimal working call needs only GroupingVar / GroupingValues. Everything else
% has a principled default. Start with defaults; tune only what diagnostics flag.

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
%     GroupingVar across all recordings of the culture. More robust for
%     multi-concentration dose-response experiments (e.g. 0, 1, 10, 100 µM).
%     Falls back to 'two_window' if fewer than Bootstrap.MinRecordings are available.
%   'metadata': read labels directly from UnitTable — see Section 6.
params.Bootstrap.GroundTruthMethod = 'two_window';

% Alpha: 1e-10 is deliberately strict. Louvain community detection handles
%   impure candidates, but over-inclusion here (weak alpha) floods the community
%   with noise units and degrades enrichment-based detection. If you get fewer
%   than ~10 responsive units per culture, consider relaxing to 1e-6 first.
%
% FullCurveAlpha: equivalent threshold for 'full_curve' (Spearman p-value).
%   0.05 is appropriate since the Spearman test is already conservative for
%   monotonic dose-response and the small N (usually 4-6 recordings).
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
% feature matrix. GroupingValues = 0 selects baseline/vehicle recordings — you
% want features from the untreated state so drug-response changes do not
% contaminate the feature space.
params.UMAP.GroupingVar    = 'Concentration';
params.UMAP.GroupingValues = 0;

% ── Unsupervised UMAP geometry ─────────────────────────────────────────────────
%
% The unsupervised UMAP is the only embedding in this pipeline. It covers all N
% unique units and its NxN fuzzy-simplicial-set graph is used for both Louvain
% community detection (generateTrainLabels) and label propagation (classifyUnits).
%
% NDims: dimensionality of the embedding. Higher values (5–10) preserve more
%   manifold structure and improve graph-based classification. 2 is interpretable
%   visually but loses information. Default 5 is a good balance.
%
% NNeighbors: controls the local/global balance of the manifold. Smaller values
%   emphasize fine-grained local structure (good for detecting small inhibitory
%   subclusters); larger values capture broader topology. AutoNNeighbors = true
%   sets this to max(MinNNeighbors, sqrt(N)) automatically — recommended.
%
% MinDist / Spread: control embedding compactness. MinDist < Spread is required.
%   These mostly affect visualization; the graph-based classifier uses the NxN
%   fuzzy simplicial set which is determined by NNeighbors, not MinDist/Spread.
%   However, BayOpt optimizes them jointly because they affect which units end up
%   in the same Louvain community.
%
% ACGWeight / WaveformWeight: L2 contribution weight of each feature group.
%   ACGs typically have ~1000 dimensions vs ~240 for waveforms; without weighting
%   ACG dominates >80% of pairwise distances. BayOpt finds the optimal balance
%   (Section 8). Default 1.0 = equal group contribution after normalization.
params.UMAP.NDims         = 5;
params.UMAP.AutoNNeighbors  = true;
params.UMAP.MinNNeighbors   = 15;
params.UMAP.MinDist         = 0.1;
params.UMAP.Spread          = 1.0;
params.UMAP.ACGWeight       = 1.0;
params.UMAP.WaveformWeight  = 1.0;

% ── Louvain community detection ────────────────────────────────────────────────
%
% Louvain community detection runs on the NxN UMAP fuzzy simplicial set.
% Communities are evaluated by their fraction of drug-responsive units to
% identify inhibitory communities.
%
% LouvainResolution (gamma): controls community granularity.
%   Higher → more, smaller communities.
%   Lower  → fewer, larger communities.
%   Default 1.0 is the standard Louvain modularity. For typical datasets
%   (100–2000 units, 5–15% responsive), the BayOpt (Section 8) will set this
%   automatically. If running without BayOpt:
%     - Raise to 1.5–3.0 if responsive units split across too many communities
%       (community plot shows many small red clusters).
%     - Lower to 0.5–0.8 if the entire dataset collapses into 2–3 communities
%       with no clear inhibitory enrichment.
%
% InhibitoryCommunityRelThresh: a community is "inhibitory" if its responsive
%   fraction >= InhibitoryCommunityRelThresh × max_community_fraction.
%   Default 0.3 means a community must have at least 30% of the highest
%   responsive fraction in the dataset. Decrease if too few communities are
%   selected; increase if CE communities are contaminated with responsive units.
%
% EnrichmentFactor: dual criterion (absolute floor). A community must have at
%   least EnrichmentFactor × (n_responsive / N_all) responsive fraction.
%   Default 3 = 3× above the dataset baseline rate. This prevents weakly-
%   enriched communities from qualifying as inhibitory when baseline p_resp is
%   high. Decrease to 2 if too few inhibitory communities are detected.
%
% PuritySigmaThreshold: responsive units with anomalously low UMAP-neighborhood
%   purity (< median - PuritySigmaThreshold × robust_sigma) are removed as
%   outliers before training label assignment. Default 2.5.
%
% CommunityFallbackThreshold: if the fraction of responsive units that land in
%   inhibitory communities is below this value, the pipeline falls back to the
%   legacy distance-based CE selection. Default 0.5 (50% coverage required).
%
% LouvainRestarts: number of Louvain runs per evaluation; the run with highest
%   modularity Q is kept. Default 5. Increase to 10 on datasets where Q is
%   unstable across runs (check Q value in diagnosticTrainLabels plot 1,2).
params.Community.LouvainResolution            = 1.0;
params.Community.InhibitoryCommunityRelThresh = 0.3;
params.Community.EnrichmentFactor            = 3;
params.Community.PuritySigmaThreshold        = 2.5;
params.Community.CommunityFallbackThreshold  = 0.5;
params.Community.LouvainRestarts             = 5;

% ── Classification ─────────────────────────────────────────────────────────────
%
% Method: "graph" (default) — label propagation on the UMAP NxN graph.
%   The graph is the same fuzzy simplicial set used for Louvain; propagation
%   is fully transductive (no projection of held-out data needed). Training
%   labels are clamped at each iteration; propagation stops when max label
%   change < GraphConvergenceTol or GraphMaxIter is reached.
%
%   "knn" — distance-weighted kNN in feature space. Simpler and faster but
%   ignores the manifold structure. Use as a sanity check or when the graph
%   is very sparse (very small datasets).
%
% AutoConfidenceK / ConfidenceK: for "knn" method only. AutoConfidenceK = true
%   sets k = max(5, sqrt(N_train)) — prevents k from approaching N_train on
%   small datasets. Ignored in "graph" mode.
%
% UseConfidenceThreshold: mark predictions below ConfidenceThreshold as NaN.
%   Default false (classify all units). Enable after inspecting the confidence
%   distribution in diagnosticClassification (panel 1,1). A threshold of 0.6
%   is a reasonable starting point for "graph" mode — values below 0.5 are
%   near-random for binary classification.
params.Classification.Method               = "graph";
params.Classification.GraphMaxIter         = 100;
params.Classification.GraphConvergenceTol  = 1e-4;
params.Classification.AutoConfidenceK      = true;    % for "knn" method
% params.Classification.UseConfidenceThreshold = true;
% params.Classification.ConfidenceThreshold    = 0.6;

% ── Counterexample selection ───────────────────────────────────────────────────
%
% In the community path (default), excitatory counterexamples are selected by
% k-medoids from non-responsive units in CE (non-inhibitory) communities.
% CounterexampleRatio = 1 gives one CE per inhibitory candidate — a balanced
% training set. Increase to 2–3 if the inhibitory fraction after classification
% is persistently > 30%; this helps the graph propagation discriminate better.
params.OutlierDetection.CounterexampleRatio = 1;

% ── Reproducibility ───────────────────────────────────────────────────────────
%
% RNGSeed controls UMAP and bootstrap random states. Run with 2–3 different seeds
% to verify stability before finalizing results. If labels flip substantially
% (ARI < 0.85), run optimizeUnsupervisedUMAP() — the embedding geometry is
% sensitive to initial conditions at the default parameter values.
params.RNGSeed = 42;

% ── Diagnostics ───────────────────────────────────────────────────────────────
%
% When enabled, each pipeline stage auto-generates a diagnostic figure:
%   identifyResponsiveUnits  → diagnosticResponsiveUnits   (responsive unit check)
%   generateTrainLabels      → diagnosticTrainLabels        (label quality check)
%   classifyUnits            → diagnosticClassification     (prediction plausibility)
%   optimizeUnsupervisedUMAP → diagnosticOptimization       (BayOpt convergence)
%
% Each can also be called manually: ctc.diagnosticTrainLabels()
%
% params.Diagnostics.Enable      = true;
% params.Diagnostics.SaveDir     = '/path/to/output';
% params.Diagnostics.ShowFigures = true;   % false suppresses display in batch mode

% Construct classifier
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

fprintf('Inhibitory candidates: %d / %d total units (%.1f%%)\n', ...
    sum(ctc.ResponsiveUnitIdx), numel(ctc.ResponsiveUnitIdx), ...
    100 * mean(ctc.ResponsiveUnitIdx));

%% 4  Generate training labels (Louvain community detection)
%
% Internally:
%   1. buildNormalizedFeatures — extracts harmonized waveforms + ACGs; builds
%      the feature matrix; applies per-group z-score (NormalizationVar groups);
%      stores NormalizationParams for classifyUnits.
%   2. Unsupervised UMAP on ALL unique units — produces an NxN fuzzy simplicial
%      set stored as ctc.UMAP.graph and a 2D (or NDims-D) embedding in
%      ctc.Reduction.Unsupervised. If ctc.UMAP is already populated (e.g. after
%      optimizeUnsupervisedUMAP), this step is skipped automatically.
%   3. Louvain community detection on ctc.UMAP.graph — identifies communities
%      by responsive-unit enrichment (dual criterion: relative + absolute).
%   4. Inhibitory training set — responsive units in inhibitory communities,
%      after graph-purity outlier removal.
%   5. Excitatory counterexamples — k-medoids selection from CE communities,
%      proportional to community size, total = n_clean_inhibitory.

ctc.generateTrainLabels();

tl = ctc.TrainLabels;
fprintf('Train set: %d excitatory, %d inhibitory (Q=%.3f)\n', ...
    sum(tl.sorted_y_train == 1), sum(tl.sorted_y_train == 2), tl.Q_modularity);

% If the community path was used, inspect community structure:
if tl.use_community_path
    fprintf('Inhibitory communities: %s\n', num2str(tl.inh_comm_ids));
end

%% 5  Classify all units (graph label propagation)
%
% Uses the NxN UMAP graph stored in ctc.UMAP.graph. Training labels are clamped
% at each iteration; label distributions propagate through edge weights until
% convergence. No test-set projection is performed — the embedding covers all
% units from the start (transductive learning).
%
% Confidence = winning class probability after propagation (0.5 = near-random,
% 1.0 = unanimous neighborhood). Training units are always confidence = 1.0.

ctc.classifyUnits();

labels = ctc.UnitLabels;   % 1=excitatory, 2=inhibitory, NaN=unclassified
n_exc  = sum(labels == 1, 'omitnan');
n_inh  = sum(labels == 2, 'omitnan');
n_nan  = sum(isnan(labels));
fprintf('Excitatory: %d  Inhibitory: %d  Unclassified: %d\n', n_exc, n_inh, n_nan);

% Typical inhibitory fraction in MEA cultures: 15–25%. Values outside this
% range suggest parameter issues — check diagnosticClassification (panel 1,2:
% per-chip I/E fraction) and diagnosticTrainLabels (panel 2,2: community fractions).

%% 6  Metadata-driven labels (pure culture scenario)
%
% When ground truth is externally supplied (optogenetics, genetic markers,
% morphology tags), use GroundTruthMethod = 'metadata' to read labels
% directly from the UnitTable. This skips the bootstrap entirely.
%
% LabelField:           column in fs.UnitTable containing cell type labels.
% ResponsiveClassValue: the string value that denotes the inhibitory class.

params_meta = params;
params_meta.Bootstrap.GroundTruthMethod    = 'metadata';
params_meta.Bootstrap.LabelField           = 'CellType';   % column in UnitTable
params_meta.Bootstrap.ResponsiveClassValue = 'inhibitory';

ctc_meta = CellTypeClassifier(fs, ud, params_meta);
ctc_meta.identifyResponsiveUnits();   % reads labels from UnitTable, returns immediately
ctc_meta.generateTrainLabels();
ctc_meta.classifyUnits();

fprintf('Metadata path: %d exc, %d inh\n', ...
    sum(ctc_meta.UnitLabels==1,'omitnan'), ...
    sum(ctc_meta.UnitLabels==2,'omitnan'));

%% 7  Inspect harmonized features
%
% HarmonizedWaveforms and HarmonizedACGs are populated by generateTrainLabels.

wf  = ctc.HarmonizedWaveforms;   % (N_samples x N_units)
acg = ctc.HarmonizedACGs;         % (N_bins   x N_units)
sr  = ctc.HarmonizedSR;

fprintf('Waveform : %d x %d at %.0f Hz\n', size(wf,1), size(wf,2), sr);
fprintf('ACG      : %d x %d\n', size(acg,1), size(acg,2));

plotCellTypeFeatures(ctc);
sortACGsByPeak(ctc);

%% 8  Bayesian optimization of UMAP + community parameters (Phase 1)
%
% Run this when default parameters produce unsatisfying community structure:
%   - diagnosticTrainLabels panel (1,2): responsive units scattered across many
%     communities with no clear inhibitory enrichment.
%   - diagnosticTrainLabels panel (2,2): community bar chart shows no community
%     with clearly elevated responsive fraction (no red bar far from the pack).
%   - Inhibitory fraction after classifyUnits is implausibly low (<10%) or high (>30%).
%   - Labels are unstable across RNG seeds (ARI < 0.85 in assessStability).
%
% What is optimized (UnsupOptimizeVars — all six by default):
%   NNeighbors       — local/global balance of the manifold. Smaller → finer
%     local communities; larger → broader topology. Interacts with LouvainResolution.
%   MinDist          — embedding compactness. Mainly visual; affects Louvain
%     indirectly through community shape.
%   Spread           — embedding spread (tune jointly with MinDist).
%   ACGWeight        — L2 weight of the ACG feature group. ACGs dominate by
%     dimension count; BayOpt finds the balance where inhibitory units form
%     a coherent community without ACG domination.
%   WaveformWeight   — L2 weight of the waveform feature group.
%   LouvainResolution — community granularity. AutoLouvainRange = true (default)
%     sets the search range to [min, N_all/(10 × n_responsive)], scaling to the
%     sparsity of responsive units — sparser datasets need finer communities.
%
% Objective: community_coherence = fraction of responsive units that land in
%   inhibitory communities (dual-criterion qualified). Higher is better.
%   Equivalently, loss = -(coherence); best_objective will be near -1.0 for
%   well-separable datasets.
%
% This is a one-time calibration step for a new experimental preparation.
% The optimized parameters are written back to ctc.Parameters automatically.

ctc_opt = CellTypeClassifier(fs, ud, params);
ctc_opt.identifyResponsiveUnits();

results = ctc_opt.optimizeUnsupervisedUMAP();
fprintf('Phase 1 BayOpt: best coherence = %.3f\n', -results.bestObjective);
fprintf('Best parameters:\n');
disp(results.bestParams);

% generateTrainLabels reuses the UMAP stored by optimizeUnsupervisedUMAP
% (ReuseExistingUMAP=true default) — no second UMAP run needed.
ctc_opt.generateTrainLabels();
ctc_opt.classifyUnits();

% Verify stability across seeds
stability = ctc_opt.assessStability('NRuns', 5);
fprintf('Stability: ARI = %.3f +/- %.3f\n', stability.meanARI, stability.stdARI);

% To fix only a subset of variables (e.g. keep NNeighbors fixed after manual tuning):
%   ctc_opt.Parameters.BayesianOptimization.UnsupOptimizeVars = ...
%       ["MinDist", "Spread", "ACGWeight", "WaveformWeight", "LouvainResolution"];
%   ctc_opt.optimizeUnsupervisedUMAP();

%% 9  Validate training labels
%
% Leave-one-culture-out cross-validation on the training set. Trains on all
% cultures except one, classifies held-out units, and reports accuracy.
% Low per-culture accuracy (<75%) indicates that culture's training labels
% are inconsistent with the rest — check diagnosticTrainLabels for that culture.

ctc.validateTrainingLabels();

%% 10  Inspect UMAP embedding and classification
%
% The unsupervised UMAP embedding covers all unique units and is stored in
% ctc.Reduction.Unsupervised. Colour by label or confidence to verify structure.

unsup = ctc.Reduction.Unsupervised;

% Classification map: excitatory blue, inhibitory red, NaN grey
labels_unique = ctc.UnitLabels(ctc.NormalizedFeatures.unique_to_rep);
conf_unique   = ctc.UnitConfidence(ctc.NormalizedFeatures.unique_to_rep);

exc_mask = labels_unique == 1;
inh_mask = labels_unique == 2;
unc_mask = isnan(labels_unique);

figure;
subplot(1, 2, 1);
hold on;
if any(unc_mask)
    scatter(unsup(unc_mask,1), unsup(unc_mask,2), 8, [.6 .6 .6], 'filled', ...
        'MarkerFaceAlpha', 0.4, 'DisplayName', 'Unclassified');
end
scatter(unsup(exc_mask,1), unsup(exc_mask,2), 15, [0.1 0.3 0.8], 'filled', ...
    'MarkerFaceAlpha', 0.7, 'DisplayName', 'Excitatory');
scatter(unsup(inh_mask,1), unsup(inh_mask,2), 15, [0.8 0.1 0.1], 'filled', ...
    'MarkerFaceAlpha', 0.7, 'DisplayName', 'Inhibitory');
hold off;
legend('Box', 'off');
title('Classification labels');
xlabel('UMAP 1'); ylabel('UMAP 2');

% Confidence map: marker size ∝ confidence, color = class
subplot(1, 2, 2);
msz_exc = max(5, 5 + 40 * conf_unique(exc_mask));
msz_inh = max(5, 5 + 40 * conf_unique(inh_mask));
hold on;
scatter(unsup(exc_mask,1), unsup(exc_mask,2), msz_exc(:), [0.1 0.3 0.8], 'filled', ...
    'MarkerFaceAlpha', 0.7, 'DisplayName', 'Excitatory');
scatter(unsup(inh_mask,1), unsup(inh_mask,2), msz_inh(:), [0.8 0.1 0.1], 'filled', ...
    'MarkerFaceAlpha', 0.7, 'DisplayName', 'Inhibitory');
hold off;
title('Classification confidence (size \propto conf)');
xlabel('UMAP 1'); ylabel('UMAP 2');

% Signs of good classification:
%   - Two spatially separated regions in the unsupervised UMAP.
%   - Inhibitory units cluster with the inhibitory community from generateTrainLabels.
%   - Training units (ctc.TrainLabels.sorted_train_ids) sit inside the relevant cluster.
%   - High confidence (large markers) in cluster cores; lower confidence at boundaries.
%
% If inhibitory units are spatially intermixed with excitatory units:
%   - Inspect diagnosticTrainLabels panel (2,2): are there clearly enriched communities?
%   - If not, run optimizeUnsupervisedUMAP() (Section 8).
%   - If yes but classification still fails, lower Community.EnrichmentFactor
%     or adjust Community.InhibitoryCommunityRelThresh.

%% 11  Recommended troubleshooting workflows
%
% ── Problem: fewer than ~10 responsive units per culture ───────────────────────
%   1. Relax Bootstrap.Alpha to 1e-6 (keep Direction = 'increase').
%   2. If still too few, switch to Bootstrap.GroundTruthMethod = 'full_curve'
%      for dose-response designs.
%   3. As last resort, use 'metadata' if external labels are available.
%
% ── Problem: inhibitory fraction < 10% or > 30% ───────────────────────────────
%   1. Check diagnosticTrainLabels panel (2,2): community responsive fractions.
%      - If the bar chart shows no clearly inhibitory community (all bars similar
%        height): run optimizeUnsupervisedUMAP to improve community structure.
%      - If an inhibitory community exists but fraction is still wrong: inspect
%        diagnosticClassification panel (1,2) per-chip I/E bars. Outlier chips
%        may indicate poor recording quality for those wells.
%   2. Increase Bootstrap.Alpha (stricter) to reduce contamination in the
%      responsive set if the fraction is too high.
%   3. Increase CounterexampleRatio (e.g. to 2) if the fraction is too high —
%      a larger excitatory training set helps the graph boundary discriminate.
%
% ── Problem: unstable results across RNG seeds ────────────────────────────────
%   1. Run ctc.assessStability('NRuns', 5). If ARI < 0.85, the embedding is
%      sensitive to initialization — run optimizeUnsupervisedUMAP.
%   2. After BayOpt, re-run assessStability. ARI > 0.90 is the target.
%   3. If still unstable: increase LouvainRestarts (e.g. to 10) or reduce
%      the BayOpt run by fixing WaveformWeight/ACGWeight at their best values
%      and only optimizing geometry parameters.
%
% ── Problem: Louvain fallback triggered ("Community fallback" warning) ─────────
%   The pipeline fell back to distance-based CE selection because fewer than
%   CommunityFallbackThreshold (50% by default) of responsive units landed in
%   inhibitory communities. This means the UMAP graph does not separate
%   inhibitory units into coherent communities at the current resolution.
%   1. Run optimizeUnsupervisedUMAP — community coherence is the direct objective.
%   2. If already optimized: lower Community.CommunityFallbackThreshold to 0.3
%      to be more permissive, or lower Community.EnrichmentFactor to 2.
%
% ── Problem: generateTrainLabels reruns UMAP despite prior optimizeUnsupervisedUMAP
%   Should not happen with default ReuseExistingUMAP = true. If it does:
%   - Check that ctc is the same object (not a new CellTypeClassifier instance).
%   - Call explicitly: ctc.generateTrainLabels('ReuseExistingUMAP', true).
%   - If you changed UMAP parameters after BayOpt: ReuseExistingUMAP=false is
%     correct behaviour — you want a fresh UMAP reflecting the new parameters.
