% run_clf_metadata.m
%
% Cell-type classification using metadata-defined training labels.
%
% Use when ground-truth or proxy labels are available in recording metadata
% (e.g. genetic markers, sorted populations, manually curated labels).
% Bypasses the bootstrap firing-rate test and unsupervised UMAP label generation.

%% 0 — Setup

deephys_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(deephys_root));

%% 1 — Configuration

proc_paths = {
    '/path/to/proc1.mat'
    '/path/to/proc2.mat'
};
save_dir = '/path/to/save';

% Metadata field that carries cell-type labels
label_field = 'CellType';   % e.g. 'CellType', 'GeneticLabel'
inh_value   = 'PV';         % value indicating inhibitory units
exc_value   = 'PYR';        % value indicating excitatory units
balance     = true;          % balance class sizes by subsampling

%% 2 — Load data

procs = RecordingProcessor.loadMany(proc_paths);
fs    = FeatureStore.fromProcessors(procs);

ud = [];
for i = 1:numel(procs)
    ud = [ud, procs(i).Units]; %#ok<AGROW>
end

%% 3 — Initialise CellTypeClassifier

parameters = struct();
parameters.UMAP.NormalizationVar   = 'ChipID';
parameters.UMAP.GroupingVar        = 'Concentration';
parameters.UMAP.GroupingValues     = 0;
parameters.UMAP.NNeighbors         = 100;
parameters.UMAP.MinDist            = 0.1;
parameters.UMAP.NDims              = 10;
parameters.Harmonization.ACGSource = 'FullACG';

ctc = CellTypeClassifier(fs, ud, parameters);

%% 4 — Build training labels from metadata

% Look up label_field in UnitTable for each unit
if ~ismember(label_field, string(fs.UnitTable.Properties.VariableNames))
    error('Field "%s" not found in UnitTable. Available: %s', ...
        label_field, strjoin(string(fs.UnitTable.Properties.VariableNames), ', '));
end

meta_vals = string(fs.UnitTable.(label_field));
inh_mask  = meta_vals == string(inh_value);
exc_mask  = meta_vals == string(exc_value);

fprintf('Inhibitory candidates: %d\n', sum(inh_mask));
fprintf('Excitatory candidates: %d\n', sum(exc_mask));
assert(sum(inh_mask) > 0 && sum(exc_mask) > 0, ...
    'No training units found — check label_field / inh_value / exc_value.');

inh_idx = find(inh_mask);
exc_idx = find(exc_mask);

if balance
    n_min   = min(numel(inh_idx), numel(exc_idx));
    rng(42);
    inh_idx = randsample(inh_idx, n_min, false);
    exc_idx = randsample(exc_idx, n_min, false);
    fprintf('Balanced: %d inhibitory + %d excitatory\n', n_min, n_min);
end

train_ids = [exc_idx(:)', inh_idx(:)'];
y_train   = [ones(1, numel(exc_idx)), 2*ones(1, numel(inh_idx))];
[sorted_train_ids, sort_idx] = sort(train_ids, 'ascend');

labels.sorted_train_ids = sorted_train_ids;
labels.sorted_y_train   = y_train(sort_idx);
labels.umap_train_idx   = false(1, height(fs.UnitTable));
labels.umap_train_idx(sorted_train_ids) = true;
labels.umap_test_idx    = ~labels.umap_train_idx;

ctc.TrainLabels = labels;

%% 5 — Classify all units

ctc.classifyUnits();
fprintf('Result: %d excitatory | %d inhibitory | %d unclassified\n', ...
    sum(ctc.UnitLabels == 1), sum(ctc.UnitLabels == 2), sum(isnan(ctc.UnitLabels)));

%% 6 — Inspect quality

plotCellTypeFeatures(ctc);

chip_col = string(fs.UnitTable.ChipID);
chips    = unique(chip_col, 'stable');
for c = 1:numel(chips)
    mask = chip_col == chips(c);
    ie   = sum(ctc.UnitLabels(mask) == 2) / sum(mask);
    fprintf('  ChipID %s: I/E fraction = %.2f\n', chips(c), ie);
end

%% 7 — E/I analysis

eia_params = struct();
eia_params.Activity.BinSize              = 0.01;
eia_params.Activity.SecCutout            = [0, 1200];
eia_params.BurstDetection.Threshold      = 2;
eia_params.BurstDetection.SmoothWindow   = 9;
eia_params.BurstDetection.PeakCutout     = 150;
eia_params.BurstDetection.PeakHeight     = 0.1;
eia_params.BurstDetection.PeakProminence = 0.1;
eia_params.BurstDetection.MinPeakDistance = 6;
eia_params.BurstDetection.SelectedVariant = 4;

eia = EIAnalyzer(ctc, eia_params);
eia.computeActivity();
eia.detectBursts();
eia.extractBurstCutouts();
eia.computeCorrelations();
eia.normalizeBurstCutouts();

eia.PlotNetworkActivity(1);

norm_bursts = eia.NormalizedCutouts.bursts;
norm_ie     = eia.NormalizedCutouts.ie;

figure('Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact');
nexttile; imagesc(norm_ie);     title('I/E ratio (normalised)',      'FontWeight', 'normal'); colorbar
nexttile; imagesc(norm_bursts); title('Total activity (normalised)', 'FontWeight', 'normal'); colorbar
nexttile; plot(nanmean(norm_ie), 'r'); hold on; plot(nanmean(norm_bursts), 'k');
          legend('I/E', 'Total'); xlabel('Bin'); box off
nexttile; plot(gradient(nanmean(norm_bursts)));
          xlabel('Bin'); title('d/dt total activity', 'FontWeight', 'normal'); box off

%% 8 — Save

if ~isfolder(save_dir); mkdir(save_dir); end
save_path = fullfile(save_dir, ...
    strjoin(["clf_metadata", label_field, string(inh_value), "vs", string(exc_value)], "_"));
save(save_path, 'ctc', 'eia', 'norm_bursts', 'norm_ie', ...
     'label_field', 'inh_value', 'exc_value', 'balance');
fprintf('Saved to %s.mat\n', save_path);
