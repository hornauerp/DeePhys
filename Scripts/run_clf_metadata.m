% run_clf_metadata.m
%
% General-purpose cell-type classification using metadata-defined training labels.
%
% Unlike run_clf_transfer.m (which identifies inhibitory candidates via a
% bootstrap firing-rate response test), this script assigns training labels
% directly from recording metadata fields.  Use this when ground-truth or
% proxy cell-type information is available in the metadata — e.g. genetic
% markers, sorted cell-type fractions, or manually curated labels.
%
% Pipeline:
%   1.  Load MEArecordings and form a RecordingGroup
%   2.  Initialise CellTypeClassifier
%   3.  Build training labels from metadata
%       a. Collect all units from all cultures
%       b. Flag inhibitory / excitatory training units via metadata filters
%       c. Optionally balance classes by subsampling the larger one
%       d. Set ctc.TrainLabels directly (bypasses bootstrap + UMAP label generation)
%   4.  Classify all units via supervised UMAP
%   5.  Inspect classification quality (ACG heatmaps, I/E ratio per chip)
%   6.  E/I network analysis (EIAnalyzer)
%   7.  Visualise burst dynamics
%   8.  Save
%
% Requirements:
%   - DeePhys on path
%   - run_umap toolbox on path
%   - Chemogenetics/Helpers on path (load_cortical_cultures)

%% 0 — Setup

deephys_root = "/home/phornauer/Git/DeePhysDev/DeePhys";
addpath(genpath(deephys_root));
addpath("/home/phornauer/Git/Chemogenetics/Helpers");

%% 1 — Load recordings and form RecordingGroup

batch_ids = [1, 1, 2, 2];
week_ids  = [2, 3, 1, 2];
full_rec_array = load_cortical_cultures(batch_ids, week_ids);

rg_params.Selection.Inclusion = {{'Region', 'Cortex'}};
rg_params.Selection.Exclusion = {{'Concentration', inf}};
rg = RecordingGroup(full_rec_array, rg_params);
rg.Cultures = rg.groupCultures('RecordingDate');

%% 2 — Initialise CellTypeClassifier

parameters = struct();

% UMAP embedding — used by classifyUnits to build the feature matrix
% Set GroupingVar / GroupingValues to select which recording condition
% to use for unit features.  Use GroupingValues = [] to include all.
parameters.UMAP.UnitFeatures    = ["FullACG", "ReferenceWaveform"];
parameters.UMAP.NormalizationVar= "ChipID";
parameters.UMAP.GroupingVar     = "Concentration";
parameters.UMAP.GroupingValues  = 0;       % 0 = baseline only
parameters.UMAP.NNeighbors      = 100;
parameters.UMAP.MinDist         = 0.1;
parameters.UMAP.Spread          = 1;
parameters.UMAP.NDims           = 10;

ctc = CellTypeClassifier(rg, parameters);

%% 3 — Build training labels from metadata

% --- 3a: Specify metadata filters ---
%
% label_field   : metadata field that carries cell-type information
% inh_value     : value(s) of label_field that indicate inhibitory units
% exc_value     : value(s) of label_field that indicate excitatory units
% balance       : if true, subsample the larger class to match the smaller
%
% Example: units whose recording has Metadata.CellType == 'PV' are
% treated as inhibitory; 'PYR' units as excitatory.
%
label_field = 'CellType';     % <-- adjust to your metadata field
inh_value   = 'PV';           % <-- value(s) marking inhibitory units
exc_value   = 'PYR';          % <-- value(s) marking excitatory units
balance     = true;

% --- 3b: Collect all units and extract the label metadata field ---
unit_list = Unit.empty;
for c = 1:length(rg.Cultures)
    unit_list = [unit_list, rg.Cultures(c).Units]; %#ok<AGROW>
end
ctc.UnitList = unit_list;

meta_vals = arrayfun(@(u) u.MEArecording.Metadata.(label_field), ...
    ctc.UnitList, 'UniformOutput', false);

inh_mask = cellfun(@(v) isequal(v, inh_value), meta_vals);
exc_mask = cellfun(@(v) isequal(v, exc_value), meta_vals);

fprintf('Inhibitory training candidates: %i\n', sum(inh_mask));
fprintf('Excitatory training candidates: %i\n', sum(exc_mask));
assert(sum(inh_mask) > 0 && sum(exc_mask) > 0, ...
    'No training units found — check label_field / inh_value / exc_value.');

% --- 3c: Optional class balancing ---
inh_idx = find(inh_mask);
exc_idx = find(exc_mask);

if balance
    n_min = min(length(inh_idx), length(exc_idx));
    inh_idx = randsample(inh_idx, n_min, false);
    exc_idx = randsample(exc_idx, n_min, false);
    fprintf('Balanced training set: %i inhibitory + %i excitatory\n', n_min, n_min);
else
    fprintf('Training set (unbalanced): %i inhibitory + %i excitatory\n', ...
        length(inh_idx), length(exc_idx));
end

% --- 3d: Build TrainLabels struct (same format as generateTrainLabels) ---
train_ids = [exc_idx(:)', inh_idx(:)'];
y_train   = [ones(1, length(exc_idx)), 2*ones(1, length(inh_idx))];

[sorted_train_ids, sort_idx] = sort(train_ids, 'ascend');

labels.sorted_train_ids = sorted_train_ids;
labels.sorted_y_train   = y_train(sort_idx);
labels.umap_train_idx   = false(1, length(ctc.UnitList));
labels.umap_train_idx(sorted_train_ids) = true;
labels.umap_test_idx    = ~labels.umap_train_idx;

ctc.TrainLabels = labels;

%% 4 — Classify all units via supervised UMAP
%
% Trains supervised UMAP on the metadata-labelled training set, then
% projects all remaining units and assigns labels by nearest supervisor.

ctc.classifyUnits();

n_exc = sum(ctc.UnitLabels == 1);
n_inh = sum(ctc.UnitLabels == 2);
n_nan = sum(isnan(ctc.UnitLabels));
fprintf('Classification: %i excitatory | %i inhibitory | %i unclassified\n', ...
    n_exc, n_inh, n_nan);

%% 5 — UMAP sanity check  (optional — comment out to skip)
%
% Shows supervised train/test embeddings coloured by label.
% No unsupervised tile since generateTrainLabels was bypassed.

plotUMAPSanityCheck(ctc);

%% 6 — Inspect classification quality

% ACG heatmaps
all_acgs  = horzcat(ctc.UnitList.FullACG)';
norm_acgs = all_acgs ./ max(all_acgs, [], 2);

inh_acgs = norm_acgs(ctc.UnitLabels == 2, :);
exc_acgs = norm_acgs(ctc.UnitLabels == 1, :);

[sorted_inh_acgs, ~] = sortACGsByPeak(inh_acgs, "max");
[sorted_exc_acgs, ~] = sortACGsByPeak(exc_acgs, "max");

figure('Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; imagesc(sorted_inh_acgs); title('Inhibitory ACGs', 'FontWeight', 'normal'); colorbar
nexttile; imagesc(sorted_exc_acgs); title('Excitatory ACGs',  'FontWeight', 'normal'); colorbar

% Waveforms and ACGs split by cell type and training / test set
plotCellTypeFeatures(ctc.UnitList, ctc.UnitLabels, ctc.TrainLabels);

% I/E fraction per chip
[chip_idx, ~] = rg.combineMetadataIndices(ctc.UnitList, "ChipID");
in_count  = histcounts(chip_idx(ctc.UnitLabels == 2));
all_count = histcounts(chip_idx);
ie_per_chip = in_count ./ all_count;
fprintf('Inhibitory fraction per chip: ');
fprintf('%.2f ', ie_per_chip);
fprintf('\n');

%% 7 — Validation against ground truth  (optional — comment out to skip)
%
% Compare predicted test-set labels against an INDEPENDENT ground truth.
% NOTE: use a different metadata field from label_field used for training,
% otherwise you are comparing against the training signal itself.
% Units with NaN ground truth are silently excluded.

gt_field     = 'IndependentCellType'; % <-- different from label_field above
gt_exc_value = 'PYR';                 % <-- value indicating excitatory (class 1)
gt_inh_value = 'PV';                  % <-- value indicating inhibitory (class 2)

gt_meta   = arrayfun(@(u) u.MEArecording.Metadata.(gt_field), ...
    ctc.UnitList, 'UniformOutput', false);
gt_labels = nan(1, length(ctc.UnitList));
gt_labels(cellfun(@(v) isequal(v, gt_exc_value), gt_meta)) = 1;
gt_labels(cellfun(@(v) isequal(v, gt_inh_value), gt_meta)) = 2;

% Evaluate on test units only (training labels are known by construction)
test_mask = ctc.TrainLabels.umap_test_idx;
metrics   = validateClassification(ctc.UnitLabels(test_mask), gt_labels(test_mask));

%% 8 — E/I network analysis

eia_params = struct();
eia_params.Activity.BinSize               = 0.01;
eia_params.Activity.SecCutout             = [0, 1200];
eia_params.BurstDetection.Threshold       = 2;
eia_params.BurstDetection.SmoothWindow    = 9;
eia_params.BurstDetection.PeakCutout      = 150;
eia_params.BurstDetection.PeakHeight      = 0.1;
eia_params.BurstDetection.PeakProminence  = 0.1;
eia_params.BurstDetection.MinPeakDistance = 6;
eia_params.BurstDetection.SelectedVariant = 4;

eia = EIAnalyzer(ctc, eia_params);
eia.computeActivity();
eia.detectBursts();
eia.extractBurstCutouts();
eia.computeCorrelations();

%% 9 — Visualise

eia.PlotNetworkActivity(1);

eia.normalizeBurstCutouts();
norm_bursts = eia.NormalizedCutouts.bursts;
norm_ie     = eia.NormalizedCutouts.ie;

figure('Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact');
nexttile; imagesc(norm_ie);     title('I/E ratio (normalised)',      'FontWeight', 'normal'); colorbar
nexttile; imagesc(norm_bursts); title('Total activity (normalised)', 'FontWeight', 'normal'); colorbar
nexttile; plot(nanmean(norm_ie), 'r'); hold on; plot(nanmean(norm_bursts), 'k');
          legend('I/E', 'Total', 'Location', 'northwest'); xlabel('Bin'); box off
nexttile; plot(gradient(nanmean(norm_bursts)));
          xlabel('Bin'); title('d/dt total activity', 'FontWeight', 'normal'); box off

%% 10 — Save

save_dir = fullfile(deephys_root, 'Data', 'Classification');
if ~isfolder(save_dir); mkdir(save_dir); end
save_path = fullfile(save_dir, ...
    strjoin(["clf_metadata", label_field, string(inh_value), "vs", string(exc_value)], "_"));
save(save_path, 'ctc', 'eia', 'norm_bursts', 'norm_ie', ...
     'label_field', 'inh_value', 'exc_value', 'balance');
fprintf('Saved to %s.mat\n', save_path);
