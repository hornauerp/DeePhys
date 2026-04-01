% run_cell_type_classification.m
%
% Full cell-type classification and E/I analysis pipeline.
%
% Pipeline:
%   1. Load processors → assemble FeatureStore
%   2. Initialise CellTypeClassifier
%   3. Bootstrap firing-rate test (inhibitory candidates)
%   4. Generate training labels (UMAP + isolation forest)
%   5. Classify all units (supervised UMAP)
%   6. E/I analysis (EIAnalyzer)
%   7. Save results

%% 0 — Setup

deephys_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(deephys_root));

%% 1 — Load data

% Paths to saved RecordingProcessor .mat files (one per recording)
proc_paths = {
    '/path/to/proc1.mat'
    '/path/to/proc2.mat'
};
save_dir = '/path/to/save';

procs = RecordingProcessor.loadMany(proc_paths);
fs    = FeatureStore.fromProcessors(procs);

% Build UnitData array in FeatureStore row order
ud = [];
for i = 1:numel(procs)
    ud = [ud, procs(i).Units]; %#ok<AGROW>
end

fprintf('FeatureStore: %d units across %d recordings\n', ...
    height(fs.UnitTable), height(fs.RecordingTable));

%% 2 — Initialise CellTypeClassifier

parameters = struct();

% Bootstrap firing-rate test
parameters.Bootstrap.PreCutout  = [0, 1200];
parameters.Bootstrap.PostCutout = [6000, 7200];
parameters.Bootstrap.BinSize    = 20;
parameters.Bootstrap.NIter      = 1000;
parameters.Bootstrap.Alpha      = 1e-10;

% UMAP embedding
parameters.UMAP.NormalizationVar   = 'ChipID';
parameters.UMAP.GroupingVar        = 'Concentration';
parameters.UMAP.GroupingValues     = 0;
parameters.UMAP.NNeighbors         = 100;
parameters.UMAP.MinDist            = 0.1;
parameters.UMAP.NDims              = 10;

% ACG source: prefer full-recording (Parent_ACG*) from FeatureStore
parameters.Harmonization.ACGSource  = 'FullACG';
parameters.Harmonization.ACGBinSize = 0.0005;
parameters.Harmonization.ACGLag     = 0.1;

% Isolation forest
parameters.OutlierDetection.ContaminationFraction = 0.5;
parameters.OutlierDetection.NObsPerLearner        = 50;
parameters.OutlierDetection.DistancePercentile    = 80;

ctc = CellTypeClassifier(fs, ud, parameters);

%% 3 — Identify inhibitory candidates

% metadata_filter restricts which cultures contribute candidates
% e.g. only cultures where AAV == 128 (DREADD-inhibitory)
ctc.identifyResponsiveUnits({'AAV', 128});
fprintf('Responsive units: %d / %d\n', sum(ctc.ResponsiveUnitIdx), numel(ctc.ResponsiveUnitIdx));

%% 4 — Generate training labels

ctc.generateTrainLabels();
fprintf('Training set: %d excitatory + %d inhibitory\n', ...
    sum(ctc.TrainLabels.sorted_y_train == 1), ...
    sum(ctc.TrainLabels.sorted_y_train == 2));

%% 5 — Classify all units

ctc.classifyUnits();
fprintf('Result: %d excitatory | %d inhibitory | %d unclassified\n', ...
    sum(ctc.UnitLabels == 1), sum(ctc.UnitLabels == 2), sum(isnan(ctc.UnitLabels)));

%% 6 — Inspect classification quality

plotCellTypeFeatures(ctc);
sortACGsByPeak(ctc);

% I/E fraction per chip
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

%% 8 — Visualise burst dynamics

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

%% 9 — Save

if ~isfolder(save_dir); mkdir(save_dir); end
save_path = fullfile(save_dir, 'cell_type_classification_results');
save(save_path, 'ctc', 'eia');
fprintf('Saved to: %s.mat\n', save_path);
