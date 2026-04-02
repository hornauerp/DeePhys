% run_clf_transfer.m
%
% Cell-type classification with train/test culture split.
%
% A subset of cultures (training_culture_idx) seeds the supervised UMAP;
% the classifier is then applied to the remaining cultures.
% Results are visualised as split heatmaps and mean traces.
%
% For cross-dataset transfer (different FeatureStore / different lab),
% use ctc.applyTo(fs_target, ud_target) instead of classifyUnits — see
% Tutorial 6 (TransferLearning.m) for that workflow.

%% 0 — Setup

script_path = mfilename('fullpath');
if isempty(script_path); script_path = which(mfilename); end
deephys_root = fileparts(fileparts(script_path));
addpath(genpath(deephys_root));

%% 1 — Load data  (checkpoint-guarded — skipped on re-runs if ctc exists)

if ~exist('ctc', 'var')

    proc_paths = {
        '/path/to/proc1.mat'
        '/path/to/proc2.mat'
    };
    save_dir = '/path/to/save';

    procs = RecordingProcessor.loadMany(proc_paths);
    fs    = FeatureStore.fromProcessors(procs);

    ud = [];
    for i = 1:numel(procs)
        ud = [ud, procs(i).Units]; %#ok<AGROW>
    end

    % CellTypeClassifier parameters
    parameters = struct();
    parameters.Bootstrap.PreCutout  = [0, 1200];
    parameters.Bootstrap.PostCutout = [6000, 7200];
    parameters.Bootstrap.BinSize    = 20;
    parameters.Bootstrap.NIter      = 1000;
    parameters.Bootstrap.Alpha      = 1e-10;
    parameters.Bootstrap.Direction  = 'increase';
    parameters.UMAP.NormalizationVar   = 'ChipID';
    parameters.UMAP.GroupingVar        = 'Concentration';
    parameters.UMAP.GroupingValues     = 0;
    parameters.UMAP.NNeighbors         = 100;
    parameters.UMAP.MinDist            = 0.1;
    parameters.UMAP.NDims              = 10;
    parameters.UMAP.TargetWeight       = 0.5;
    parameters.Harmonization.ACGSource  = 'FullACG';
    parameters.Harmonization.ACGBinSize = 0.0005;
    parameters.Harmonization.ACGLag     = 0.1;
    parameters.OutlierDetection.OutlierAlpha            = 0.01;
    parameters.OutlierDetection.DipTestAlpha            = 0.05;
    parameters.OutlierDetection.MaxResponsiveComponents = 3;
    parameters.OutlierDetection.CounterexampleRatio     = 2;

    % Optional adaptive parameters (uncomment to enable):
    %   parameters.UMAP.AutoNNeighbors              = true;
    %   parameters.UMAP.AutoNDims                   = true;
    %   parameters.UMAP.AutoTargetWeight            = true;
    %   parameters.UMAP.FeatureSelection            = true;
    %   parameters.OutlierDetection.AutoCounterexampleRatio = true;

    ctc = CellTypeClassifier(fs, ud, parameters);
    ctc.identifyResponsiveUnits({'AAV', 128});
    ctc.generateTrainLabels();
    ctc.classifyUnits();

    fprintf('Classification: %d excitatory | %d inhibitory | %d unclassified\n', ...
        sum(ctc.UnitLabels == 1), sum(ctc.UnitLabels == 2), sum(isnan(ctc.UnitLabels)));

end

%% 2 — Train / test culture split

% Identify unique cultures in FeatureStore order
culture_ids    = FeatureStore.getCultureIDsForUnits( ...
    fs.UnitTable.RecordingID, fs.MetadataTable, ctc.Parameters.CultureKeys);
unique_cultures = unique(culture_ids, 'stable');
n_cultures      = numel(unique_cultures);

% Specify which culture indices are the training set
training_culture_idx = 4;   % 1-based index into unique_cultures
test_culture_idx     = setdiff(1:n_cultures, training_culture_idx);

fprintf('Training cultures: %s  |  Test cultures: %s\n', ...
    mat2str(training_culture_idx), mat2str(test_culture_idx));

%% 3 — I/E fraction per chip

chip_col = string(fs.UnitTable.ChipID);
chips    = unique(chip_col, 'stable');
for c = 1:numel(chips)
    mask = chip_col == chips(c);
    ie   = sum(ctc.UnitLabels(mask) == 2) / sum(mask);
    fprintf('  ChipID %s: I/E fraction = %.2f\n', chips(c), ie);
end

%% 4 — ACG inspection

acg_mat  = ctc.HarmonizedACGs';   % (N_units x N_bins), already normalized
inh_acgs = acg_mat(ctc.UnitLabels == 2, :);
exc_acgs = acg_mat(ctc.UnitLabels == 1, :);

figure('Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; imagesc(sortrows(inh_acgs)); title('Inhibitory ACGs', 'FontWeight', 'normal'); colorbar
nexttile; imagesc(sortrows(exc_acgs)); title('Excitatory ACGs',  'FontWeight', 'normal'); colorbar

plotCellTypeFeatures(ctc);

%% 5 — E/I analysis

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

norm_bursts = eia.NormalizedCutouts.bursts;
norm_ie     = eia.NormalizedCutouts.ie;

%% 6 — Visualise (train / test split)

sorted_idx = [training_culture_idx(:)', test_culture_idx(:)'];

figure('Color', 'w');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; imagesc(norm_ie(sorted_idx, :));     colorbar; axis tight
title('I/E ratio (train | test)', 'FontWeight', 'normal')

nexttile; imagesc(norm_bursts(sorted_idx, :)); colorbar; axis tight
title('Total activity (train | test)', 'FontWeight', 'normal')

nexttile;
plot(nanmean(norm_ie), 'r'); hold on; plot(nanmean(norm_bursts), 'k');
legend('I/E', 'Total', 'Location', 'northwest'); xlabel('Bin'); box off
title('All cultures', 'FontWeight', 'normal')

nexttile;
plot(nanmean(norm_ie(training_culture_idx, :)), 'r'); hold on
plot(nanmean(norm_bursts(training_culture_idx, :)), 'k');
title('Training', 'FontWeight', 'normal'); box off

nexttile;
plot(nanmean(norm_ie(test_culture_idx, :)), 'r'); hold on
plot(nanmean(norm_bursts(test_culture_idx, :)), 'k');
title('Test', 'FontWeight', 'normal'); box off

%% 7 — Save

if ~isfolder(save_dir); mkdir(save_dir); end
save_path = fullfile(save_dir, ...
    strjoin(["clf_transfer", string(training_culture_idx)], "_"));
save(save_path, 'ctc', 'eia', 'norm_bursts', 'norm_ie', ...
     'training_culture_idx', 'test_culture_idx');
fprintf('Saved to %s.mat\n', save_path);
