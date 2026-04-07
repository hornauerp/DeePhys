%% Tutorial 4 — Cell Type Classification
%
% Identifies excitatory (1) and inhibitory (2) neurons using a supervised
% UMAP pipeline: bootstrap firing-rate test → label generation → classification.
%
% Prerequisites:
%   - Tutorial 1 completed (FeatureStore and RecordingProcessors saved)
%   - Parent ACGs available in FeatureStore (optional but recommended)
%   - DeePhys on the MATLAB path

%% 1  Load data

fs_file   = '/path/to/FeatureStore.mat';
proc_paths = {'/path/to/proc1.mat', '/path/to/proc2.mat'};

fs    = FeatureStore.load(fs_file);
procs = RecordingProcessor.loadMany(proc_paths);

% Build UnitData array (must match UnitTable row order)
ud = [];
for i = 1:numel(procs)
    ud = [ud, procs(i).Units]; %#ok<AGROW>
end

%% 2  Construct CellTypeClassifier

% Default parameters — uses FullACG (parent recording) if available
ctc = CellTypeClassifier(fs, ud);

% With parameter overrides:
%   params.Harmonization.ACGSource = 'FullACG';   % prefer Parent_ACG* from FeatureStore
%   params.Harmonization.ACGBinSize = 0.0005;     % 0.5 ms bins
%   params.Harmonization.ACGLag     = 0.1;        % +-100 ms lag
%   params.UMAP.NormalizationVar    = 'ChipID';   % per-chip z-score before UMAP
%   ctc = CellTypeClassifier(fs, ud, params);

disp(ctc.Parameters.Harmonization);

%% 3  Bootstrap firing-rate test (inhibitory candidate identification)
%
% Compares pre-stimulus vs post-stimulus spike rates across cultures using a
% permutation test. Units with a significant rate increase are flagged as
% potential inhibitory units and used as the positive class in UMAP training.
%
% Optional: restrict which cultures contribute candidates via a metadata filter.

ctc.identifyResponsiveUnits();

% With filter (only use control cultures, concentration = 0):
%   ctc.identifyResponsiveUnits({'Concentration', 0});

fprintf('Inhibitory candidates: %d / %d\n', ...
    sum(ctc.ResponsiveUnitIdx), numel(ctc.ResponsiveUnitIdx));

%% 4  Generate training labels
%
% Runs unsupervised UMAP to embed all units, then applies adaptive outlier
% detection (dip test + Mahalanobis/GMM) on responsive candidates to form the
% inhibitory training set. Counterexamples are selected via farthest-point
% sampling from the non-responsive pool. Normalization parameters are stored
% for consistent application in classifyUnits and applyTo.

ctc.generateTrainLabels();

% Inspect training labels
tl = ctc.TrainLabels;
fprintf('Train set: %d excitatory, %d inhibitory\n', ...
    sum(tl.sorted_y_train == 1), sum(tl.sorted_y_train == 2));

%% 5  Classify all units (supervised UMAP)
%
% Trains a supervised UMAP on the labelled training set, projects all
% remaining units into the embedding, and assigns labels by nearest
% neighbour lookup.

ctc.classifyUnits();

% Results
labels = ctc.UnitLabels;   % 1=excitatory, 2=inhibitory, NaN=unclassified
n_exc  = sum(labels == 1, 'omitnan');
n_inh  = sum(labels == 2, 'omitnan');
n_nan  = sum(isnan(labels));
fprintf('Excitatory: %d  Inhibitory: %d  Unclassified: %d\n', n_exc, n_inh, n_nan);

%% 6  Inspect harmonized features
%
% After classification, harmonized waveforms and ACGs (resampled to a common
% grid) are stored for inspection and plotting.

wf  = ctc.HarmonizedWaveforms;   % (N_samples x N_units)
acg = ctc.HarmonizedACGs;         % (N_bins   x N_units)
sr  = ctc.HarmonizedSR;           % target waveform sampling rate

fprintf('Waveform size : %d x %d at %.0f Hz\n', size(wf,1), size(wf,2), sr);
fprintf('ACG size      : %d x %d\n', size(acg,1), size(acg,2));

%% 7  Plot cell type features

plotCellTypeFeatures(ctc);

%% 8  Sort ACGs by peak lag

sortACGsByPeak(ctc);

%% 9  Parameter tuning — ACG bin size and lag
%
% Changing ACGBinSize or ACGLag clears the in-memory and disk cache
% so features are recomputed on the next call to classifyUnits.

ctc2 = CellTypeClassifier(fs, ud);
ctc2.Parameters.Harmonization.ACGBinSize = 0.001;   % 1 ms bins
ctc2.Parameters.Harmonization.ACGLag     = 0.05;    % +-50 ms
ctc2.identifyResponsiveUnits();
ctc2.generateTrainLabels();
ctc2.classifyUnits();

%% 10  Adaptive parameters — data-driven tuning
%
% Instead of manually setting UMAP and outlier detection parameters, enable
% automatic data-driven selection. Each Auto* flag replaces a fixed constant
% with a value computed from the dataset itself.
%
% All Auto* flags default to false, so the standard pipeline is unchanged
% unless you opt in.

params_auto = struct();

% Adaptive UMAP neighborhood sizes — set to max(MinNNeighbors, sqrt(N))
% Prevents n_neighbors from being larger than the dataset (small experiments)
% or too local (large experiments).
params_auto.UMAP.AutoNNeighbors = true;
params_auto.UMAP.MinNNeighbors  = 15;

% Adaptive kNN confidence — set to max(5, sqrt(N_train))
% Scales confidence scoring with training set density.
params_auto.UMAP.AutoConfidenceK = true;

ctc_auto = CellTypeClassifier(fs, ud, params_auto);
ctc_auto.identifyResponsiveUnits();
ctc_auto.generateTrainLabels();
ctc_auto.classifyUnits();

fprintf('Fixed params : %d exc, %d inh\n', ...
    sum(ctc.UnitLabels==1,'omitnan'), sum(ctc.UnitLabels==2,'omitnan'));
fprintf('Auto  params : %d exc, %d inh\n', ...
    sum(ctc_auto.UnitLabels==1,'omitnan'), sum(ctc_auto.UnitLabels==2,'omitnan'));

%% 11  FDR-corrected responsive unit identification
%
% Replace the fixed alpha = 1e-10 with Benjamini-Hochberg FDR control.
% This scales correctly with the number of units being tested: large datasets
% get appropriate multiple-comparison correction while small datasets retain
% statistical power.

params_fdr = struct();
params_fdr.Bootstrap.UseFDR   = true;
params_fdr.Bootstrap.FDRLevel = 0.05;  % 5% false discovery rate

ctc_fdr = CellTypeClassifier(fs, ud, params_fdr);
ctc_fdr.identifyResponsiveUnits();

fprintf('Fixed alpha (1e-10): %d responsive\n', sum(ctc.ResponsiveUnitIdx));
fprintf('FDR (q=0.05):        %d responsive\n', sum(ctc_fdr.ResponsiveUnitIdx));

%% 12  Unsupervised feature selection
%
% Remove low-variance and highly correlated features before UMAP to improve
% embedding stability and reduce computation time. The selection mask is stored
% in NormalizationParams for consistent application in classifyUnits and applyTo.
%
% MinVariancePercentile: drop features in the bottom N% by variance.
% MaxCorrelation: drop one of each pair with |r| above this threshold.

params_feat = struct();
params_feat.UMAP.FeatureSelection      = true;
params_feat.UMAP.MinVariancePercentile = 10;   % drop bottom 10% by variance
params_feat.UMAP.MaxCorrelation        = 0.99; % drop one of each |r|>0.99 pair

ctc_feat = CellTypeClassifier(fs, ud, params_feat);
ctc_feat.identifyResponsiveUnits();
ctc_feat.generateTrainLabels();
ctc_feat.classifyUnits();

%% 13  Ensemble classification — robust labels via seed averaging
%
% Run the pipeline with multiple RNG seeds and take the majority vote.
% Eliminates sensitivity to a single seed choice.
% Confidence = fraction of seeds that agree on each unit's label.
% Units with agreement below MinAgreement become NaN (unclassified).

params_ens = struct();
params_ens.Ensemble.Enabled      = true;
params_ens.Ensemble.Seeds        = [42, 1042, 2042, 3042, 4042];
params_ens.Ensemble.MinAgreement = 0.6;  % need 60% seed agreement for label

ctc_ens = CellTypeClassifier(fs, ud, params_ens);
ctc_ens.identifyResponsiveUnits();
ctc_ens.classifyUnitsEnsemble();

fprintf('Ensemble: %d exc, %d inh, %d unclassified\n', ...
    sum(ctc_ens.UnitLabels==1,'omitnan'), ...
    sum(ctc_ens.UnitLabels==2,'omitnan'), ...
    sum(isnan(ctc_ens.UnitLabels)));

% Inspect per-unit agreement (confidence)
fprintf('Mean confidence: %.2f\n', mean(ctc_ens.UnitConfidence, 'omitnan'));

%% 14  Bayesian optimization — streamlined with adaptive parameters
%
% When Auto* flags are enabled, those parameters are excluded from the Bayesian
% search. MaxEvaluations is adjusted automatically to 10 * n_remaining_vars.
% This focuses optimization on genuinely free parameters.

params_opt = struct();
params_opt.UMAP.AutoNNeighbors  = true;
params_opt.UMAP.AutoConfidenceK = true;

ctc_opt = CellTypeClassifier(fs, ud, params_opt);
ctc_opt.identifyResponsiveUnits();

% Bayesian optimization searches over the remaining free parameters only
results = ctc_opt.optimizeHyperparameters();
fprintf('Optimized %d variables, best objective: %.4f\n', ...
    results.nVars, results.bestObjective);

ctc_opt.Parameters = parseStructParameters(ctc_opt.Parameters, results.bestParams);
ctc_opt.generateTrainLabels();
ctc_opt.classifyUnits();

% Verify stability of optimized + adaptive parameters
stability = ctc_opt.assessStability('NRuns', 5);
fprintf('Stability: ARI = %.3f +/- %.3f\n', stability.meanARI, stability.stdARI);

%% 15  Metadata-driven labels (skip bootstrap)
%
% If ground-truth labels are known (e.g. from optogenetics or marker staining),
% you can skip identifyResponsiveUnits and generateTrainLabels and assign
% TrainLabels manually.

n_units = numel(ud);
y_known = nan(1, n_units);   % fill with 1 (exc) or 2 (inh) for labelled units

% Example: first 10 units are excitatory, next 10 are inhibitory
y_known(1:10)  = 1;
y_known(11:20) = 2;

train_idx = ~isnan(y_known);

tl_manual = struct();
tl_manual.sorted_train_ids  = find(train_idx);
tl_manual.sorted_y_train    = y_known(train_idx);
tl_manual.umap_train_idx    = train_idx;
tl_manual.umap_test_idx     = ~train_idx;

ctc3 = CellTypeClassifier(fs, ud);
ctc3.TrainLabels = tl_manual;
ctc3.classifyUnits();

%% 16  Inspect UMAP embedding

figure;
subplot(1,2,1);
scatter(ctc.Reduction.Train(:,1), ctc.Reduction.Train(:,2), 10, ...
    ctc.TrainLabels.sorted_y_train, 'filled');
colormap([0.2 0.6 0.9; 0.9 0.3 0.2]);
title('Training embedding'); xlabel('UMAP 1'); ylabel('UMAP 2');

subplot(1,2,2);
test_embed = ctc.Reduction.Test;
test_label = ctc.UnitLabels(ctc.TrainLabels.umap_test_idx);
scatter(test_embed(:,1), test_embed(:,2), 10, test_label, 'filled');
colormap([0.2 0.6 0.9; 0.9 0.3 0.2]);
title('Test embedding'); xlabel('UMAP 1');
