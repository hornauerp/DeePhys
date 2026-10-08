%% Tutorial 4 — Cell Type Classification
%
% Classifies neurons as excitatory (1) or inhibitory (2) using a transductive
% graph-based pipeline. All units share one deterministic kNN graph (exact
% neighbours, UMAP's fuzzy weighting); communities that are statistically enriched
% in drug-responsive units provide the inhibitory seeds, depleted ones the
% excitatory seeds (consensus over Leiden restarts); label propagation on the
% graph assigns cell types to every unit.
%
%   Drug-response (default): units whose firing increases with an activating
%   DREADD in interneurons are the inhibitory ground-truth candidates. Features
%   are drug-free (baseline ACGs, smoothed with a width chosen without labels).
%
%   Metadata (pure cultures): Labels read directly from a UnitTable column
%   (e.g. optogenetics tag, genetic marker). Skips bootstrap entirely.
%
% Pipeline:
%   identifyResponsiveUnits()       ← flag drug-responsive (inhibitory) candidates
%   discardSegmentUnits()           ← keep baseline rows (dose-response data)
%   generateTrainLabels()           ← unit graph → enrichment seeds (inh / exc)
%   classifyUnits()                 ← label propagation on the unit graph
%   qualityReport()                 ← acceptance checks instead of parameter tuning
%
% The defaults are not meant to be tuned; check results with qualityReport
% (Section 5b) rather than adjusting parameters.
%
% Prerequisites:
%   - Tutorial 1 completed (FeatureStore and RecordingProcessors saved)
%   - DeePhys on the MATLAB path
%   - Brain Connectivity Toolbox (bundled; graph features and the "louvain" option)

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

%% 2  Parameters
%
% Nothing here needs tuning. What you set describes the experiment and the
% features; the remaining settings are fixed and the result is checked with
% qualityReport (Section 5b). CellTypeClassifier.returnDefaultParams lists which
% parameters the default pipeline reads; most others only serve the optional
% methods of Sections 8-9.

params = struct();

% ── Experiment description (must match your data) ──────────────────────────────
%
% CultureKeys: metadata columns that together identify one culture (well).
% GroupingVar / GroupingValues: dose variable and its baseline value. Features come
%   from the baseline recordings, so the drug response does not enter them.
% NormalizationVar: per-group z-score of the features. 'RecordingID' removes
%   recording-level offsets (electrodes, culture density) and suits cultures of
%   natural cell-type composition. Use '' for a dataset enriched in one class, where
%   per-group z-scoring would remove the class difference itself.
% GroundTruthMethod 'full_curve': candidates are units whose firing rate follows the
%   dose (Section 3). Direction: their response ('increase' for an activating DREADD
%   in interneurons). The cultures that provide candidates are chosen in Section 3.
% ControlFilter: non-transduced control cultures (qualityReport and
%   validateIndependentResponders; also the 'control_calibrated' responder method).
%   In this tutorial a metadata column 'Construct' codes the treatment: 0 control,
%   1 activating DREADD in interneurons, 2 silencing DREADD in interneurons.
%   Replace it with your own column and values.
params.CultureKeys                = ["ChipID", "PlatingDate", "RecordingDate"];
params.UMAP.GroupingVar            = 'Concentration';
params.UMAP.GroupingValues         = 0;
params.UMAP.NormalizationVar       = 'RecordingID';
params.Bootstrap.GroundTruthMethod = 'full_curve';
params.Bootstrap.Direction         = 'increase';
params.Bootstrap.ControlFilter     = {'Construct', 0};

% ── Feature definition (fixed by what the data contain) ────────────────────────
%
% Waveforms: trough-aligned, resampled to 120 kHz, window -1 to +2 ms
%   (WaveformPreTrough / WaveformPostTrough); shorter waveforms are edge-padded.
% ACGBinSize: 0.5 ms captures fine rebound inhibition peaks. Increase to 1 ms
%   for speed on large datasets. Do not go below 0.2 ms — bins become too sparse
%   for low-firing units.
% ACGLag: ±100 ms covers the full rebound window for most inhibitory subtypes.
% To classify external data later (Tutorial 6), train with the waveform window and
%   ACG lag that both datasets cover.
% ACGWeight / WaveformWeight (UMAP.*): equal total contribution of both feature
%   groups (1:1), independent of their number of bins.
% ACGSource: 'ACG' (default) uses the baseline recording only. With drug-response
%   labels this keeps the features drug-free: 'FullACG' (parent recording) and
%   'Concatenated' include the drug periods, so the ACGs of responsive units encode
%   the very response used to label them. Use 'FullACG' only without a perturbation.
% ACGSmoothing: 'auto' (default) chooses the Gaussian width from the data by
%   split-half prediction (no labels); the same width is applied to external units.
params.Harmonization.ACGBinSize   = 0.0005;
params.Harmonization.ACGLag       = 0.1;
params.Harmonization.ACGSource    = 'ACG';
params.Harmonization.ACGSmoothing = 'auto';

% ── Fixed settings (not tuned per dataset) ─────────────────────────────────────
%
% Responder calls ('full_curve'): |Spearman rho| >= Bootstrap.MonotonicityThreshold
%   (0.75), fold change >= MinFoldChange (1.5), bootstrap p < ConfirmationAlpha (1e-3).
% Seeds: exact kNN graph with UMAP.NNeighbors = 50; Leiden communities
%   (Community.Algorithm; "louvain" also available) at Community.LouvainResolution = 1
%   (standard modularity), LouvainRestarts = 30; Community.SeedFDR = 0.05.
% Sensitivity check: before reporting an inhibitory fraction, vary the responder
% thresholds (rho, fold change, alpha) and, separately, NNeighbors and
% LouvainResolution. Check that the labels agree across settings and that the
% independent validation AUC (Section 5b) stays stable. The ranking of the units is
% expected to be stable; the inhibitory fraction can move by a few percentage
% points, so report it with that range rather than tuning it.
%
% No longer needed (only read by the optional methods of Sections 8-9): the UMAP
% geometry (NDims, MinDist, Spread: plots only), EnrichmentFactor,
% InhibitoryCommunityRelThresh, PuritySigmaThreshold, the fallback settings,
% CounterexampleRatio, the isolation-forest settings, Ensemble and the Bayesian
% optimisation.

% ── Reproducibility ───────────────────────────────────────────────────────────
%
% The graph is exact and the seeds are a consensus over the community restarts, so
% the result should hardly depend on RNGSeed; compare the labels of two seeds to check.
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
% 'full_curve': per unit, the firing rate must follow the dose monotonically,
% change at least 1.5-fold and pass a bootstrap test of baseline vs highest dose
% (with the network-wide rate change subtracted). Units that pass become the
% inhibitory candidates. They need not be pure: Section 4 only uses graph
% communities in which candidates are over-represented.
%
% The filter argument restricts which cultures contribute candidates: only those
% where the drug acts on interneurons. In cultures with an inhibitory DREADD,
% units that increase are disinhibited excitatory neurons and must not be used.

% Candidates only from cultures where the activating DREADD targets interneurons:
ctc.identifyResponsiveUnits({'Construct', 1});
ctc.discardSegmentUnits();   % dose-response data: continue with the baseline rows

fprintf('Inhibitory candidates: %d / %d total units (%.1f%%)\n', ...
    sum(ctc.ResponsiveUnitIdx), numel(ctc.ResponsiveUnitIdx), ...
    100 * mean(ctc.ResponsiveUnitIdx));

%% 4  Generate training labels (enrichment seeds on the unit graph)
%
% Internally (defaults):
%   1. buildNormalizedFeatures — harmonized waveforms + drug-free ACGs (smoothed);
%      feature matrix; per-group z-score (NormalizationVar groups).
%   2. Unit graph — exact kNN fuzzy graph on all unique units (ctc.graphMatrix());
%      deterministic. A UMAP embedding is also computed for plots
%      (UMAP.ComputeEmbedding) but does not affect the labels.
%   3. Leiden on the graph, LouvainRestarts times; per community a hypergeometric
%      test of its responder count (BH-FDR Community.SeedFDR), among units of
%      cultures where responders occur.
%   4. Seeds by consensus over the restarts: responders in an enriched community in
%      >= half of them are inhibitory seeds; as many excitatory seeds come from
%      non-responders in depleted communities (lowest responder density).
% No EnrichmentFactor / purity / fallback thresholds; an error is raised when no
% responder is consistently enriched (the calls do not match the feature graph).

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
% Uses the same unit graph (ctc.graphMatrix()). Training labels are clamped at
% each iteration; label distributions propagate through edge weights until
% convergence. The graph covers all units from the start (transductive learning).
% A unit is inhibitory when P(inhibitory) > 0.5 after propagation; the training
% set is balanced, so 0.5 is the calibrated point within the training data.
%
% Confidence = winning class probability after propagation (0.5 = near-random,
% 1.0 = unanimous neighborhood). Training units are always confidence = 1.0.

ctc.classifyUnits();

labels = ctc.UnitLabels;   % 1=excitatory, 2=inhibitory, NaN=unclassified
n_exc  = sum(labels == 1, 'omitnan');
n_inh  = sum(labels == 2, 'omitnan');
n_nan  = sum(isnan(labels));
fprintf('Excitatory: %d  Inhibitory: %d  Unclassified: %d\n', n_exc, n_inh, n_nan);

% Judge your result with qualityReport below rather than by adjusting parameters.

%% 5b  Quality report (instead of tuning)
%
% Acceptance checks without ground truth: seed counts, responder FDR (with the
% control_calibrated method), the inhibitory fraction per treatment condition
% (cultures of one preparation should agree; a fraction that tracks the drug
% condition means the features encode the drug response), an independent
% validation, and the fraction with a 95% CI over cultures.
%
% Independent validation: cultures with an inhibitory DREADD in interneurons
% (hM4Di) were not used for training. Units that decrease under the actuator are
% silenced interneurons, units that increase are disinhibited excitatory neurons;
% the AUC of P(inhibitory) between the two groups measures the classification.
report = ctc.qualityReport(ConditionVar="Construct", ValidationFilter={'Construct', 2});
% or separately: v = ctc.validateIndependentResponders({'Construct', 2});

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
sortACGsByPeak(ctc.HarmonizedACGs');

%% 8  Optional (legacy): Bayesian optimization of graph + community parameters
%
% Not needed with the defaults: the labels do not depend on the UMAP geometry
% (MinDist, Spread, NDims; plots only), hardly depend on the RNG seed, and the
% enrichment seeds do not use the community_coherence objective optimised here.
% The optimiser writes NNeighbors, LouvainResolution and the feature weights back
% to ctc.Parameters; see the sensitivity check in Section 2.
% Kept for exploring Graph.Method "umap" or the "community" / "iforest" seeds.
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
ctc_opt.identifyResponsiveUnits({'Construct', 1});
ctc_opt.discardSegmentUnits();

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

%% 9  Optional: cross-validate the training labels
%
% Leave-one-culture-out cross-validation on the training set. Trains on all
% cultures except one, classifies held-out units, and reports accuracy.
% Note: it uses a supervised UMAP classifier, not the graph propagation of
% Section 5; the independent hM4Di check in qualityReport tests the actual labels.
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
% The embedding is for inspection only; the labels come from the unit graph.

%% 11  Troubleshooting (follow the qualityReport checks, not parameters)
%
% ── generateTrainLabels: "no community is enriched in responders" ─────────────
%   The responder calls do not match the feature graph.
%   1. diagnosticResponsiveUnits: are there candidates in the expected cultures
%      (identifyResponsiveUnits filter, Direction)?
%   2. Features: GroupingValues must select the baseline recordings, and
%      NormalizationVar the right groups (Section 2).
%   3. Few candidates per culture: GroundTruthMethod 'control_calibrated' with
%      Bootstrap.TargetFDR (e.g. 0.2) calls more responders at a stated purity.
%
% ── [WARN] Condition balance ──────────────────────────────────────────────────
%   The inhibitory fraction tracks the treatment, so the features encode the drug
%   response: keep ACGSource = 'ACG' (baseline only) and check GroupingValues.
%
% ── [WARN] Independent validation ─────────────────────────────────────────────
%   Expected interneurons (hM4Di decreasers) are not ranked above expected
%   excitatory units. Check the candidate cultures and Direction first; then
%   whether the validation cultures were recorded like the training cultures.
%
% ── Labels differ between RNG seeds ───────────────────────────────────────────
%   Should not happen with Graph.Method "exact" and "enrichment" seeds. If you use
%   "umap" or the embedding-based seed methods, enable Ensemble or optimise them
%   (Section 8).
%
% ── External data: (almost) all units in one class ────────────────────────────
%   See Tutorial 6: matching window, normalisation, and Cutoff="auto".
