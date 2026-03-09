% run_cell_type_classification.m
%
% Cell-type classification and E/I network analysis using DeePhys.
%
% This script recapitulates the workflow originally implemented in
% U:\Git\Chemogenetics\Scripts\Classification\clf_transfer.m and
% clf_all_batches.m, using the new @CellTypeClassifier and @EIAnalyzer
% classes. Each section notes the equivalent old function call.
%
% Pipeline overview:
%   1.  Load MEArecordings and form a RecordingGroup
%   2.  Initialise CellTypeClassifier with parameters
%   3.  Identify inhibitory candidates (bootstrap firing-rate test)
%   4.  Generate training labels (UMAP + isolation forest)
%   5.  Classify all units (supervised UMAP)
%   6.  Inspect classification quality (ACG plots, I/E ratio per chip)
%   7.  Initialise EIAnalyzer and run E/I analysis
%   8.  Visualise burst dynamics
%   9.  Save results
%
% Requirements:
%   - DeePhys on path (run addpath lines below)
%   - run_umap toolbox on path (MATLAB File Exchange)
%   - sorting outputs from Kilosort (spike_times.npy etc.)

%% 0 — Setup

deephys_root = "/home/phornauer/Git/DeePhysDev/DeePhys"; %fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(deephys_root));
addpath("/home/phornauer/Git/Chemogenetics/Helpers")

%% 1 — Load recordings and form RecordingGroup
batch_ids = 1; 
week_ids = 3;
% Point to your sorting output directories (one per recording segment).
full_rec_array = load_cortical_cultures(batch_ids, week_ids);
%      OR recording_array_from_single_files(sorting_path_list, rm_con)
%%

% rm_con       = false;    % set true to remove connectivity matrices to save memory
% full_rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

% Filter to the recordings of interest.
% Inclusion: must be cortical, AAV 0 (control), 128 (DREADD-excitatory), or 129 (DREADD-inhibitory).
% Exclusion: drop infinite-concentration (washout) entries.
rg_params.Selection.Inclusion = {{'Region','Cortex'}, {'AAV', 0, 128, 129}};
rg_params.Selection.Exclusion = {{'Concentration', inf}};

rg = RecordingGroup(full_rec_array, rg_params);

% Group recordings into Culture objects by recording date (= by dose-response series).
% old: dr_rg.Cultures = dr_rg.groupCultures('RecordingDate')
% new: Cultures is now a Culture array (not a cell array); rg.Cultures(c) not {c}
rg.Cultures = rg.groupCultures('RecordingDate');

%% 2 — Initialise CellTypeClassifier

% All parameters below match the values used in clf_transfer.m /
% clf_all_batches.m.  Override only what differs for your experiment.
parameters = struct();

% Bootstrap firing-rate test (Step 3)
parameters.Bootstrap.PreCutout  = [0, 1200];    % baseline window (s)
parameters.Bootstrap.PostCutout = [6000, 7200]; % post-DCZ window (s)
parameters.Bootstrap.BinSize    = 20;           % bin width (s)
parameters.Bootstrap.NIter      = 1000;         % permutation iterations
parameters.Bootstrap.Alpha      = 1e-10;        % two-tailed significance level

% UMAP embedding (Steps 4–5)
parameters.UMAP.UnitFeatures    = ["FullACG", "ReferenceWaveform"];
parameters.UMAP.NormalizationVar= "ChipID";     % normalise within chip
parameters.UMAP.GroupingVar     = "Concentration";
parameters.UMAP.GroupingValues  = 0;            % use baseline recordings only
parameters.UMAP.NNeighbors      = 100;
parameters.UMAP.MinDist         = 0.1;
parameters.UMAP.Spread          = 1;
parameters.UMAP.NDims           = 10;

% Isolation forest (Step 4)
parameters.OutlierDetection.ContaminationFraction = 0.5;
parameters.OutlierDetection.NObsPerLearner        = 50;
parameters.OutlierDetection.DistancePercentile    = 80;

ctc = CellTypeClassifier(rg, parameters);

%% 3 — Identify inhibitory candidates via bootstrap firing-rate test
%
% For each AAV=128 culture, units with a significant firing-rate increase
% between the pre-CNO baseline and the post-CNO window are flagged as
% inhibitory candidates.
%
% old: [interneuron_idx, unit_list] = identify_interneurons(dr_rg)
%      which internally called bootstrap_fr_response() per culture
metadata_filter = {'AAV', 128};
ctc.identifyResponsiveUnits(metadata_filter);

%% 4 — Generate training labels
%
% 1. Runs unsupervised UMAP on all units (FullACG + waveform, baseline only).
% 2. Applies isolation forest to the candidate cluster to remove outliers.
% 3. Selects counterexamples: units farther than the 95th-percentile
%    distance from the candidate centroid, resampled to match candidate count.
% 4. Builds a balanced training set (class 1 = excitatory, class 2 = inhibitory).
%
% old: [labels, reduction] = generate_train_labels(dr_rg, unit_list,
%                              interneuron_idx, umap_params, outlier_params)

ctc.generateTrainLabels();

n_inh_train = sum(ctc.TrainLabels.sorted_y_train == 2);
n_exc_train = sum(ctc.TrainLabels.sorted_y_train == 1);
fprintf('Training set: %i inhibitory + %i excitatory examples\n', ...
    n_inh_train, n_exc_train);

%% 5 — Classify all units via supervised UMAP
%
% Trains a supervised UMAP on the labelled training set, then projects
% held-out units and assigns labels by nearest supervisor.
%
% old: full_umap_labels = classify_units(dr_rg, unit_list, labels, umap_params)

ctc.classifyUnits();

n_exc = sum(ctc.UnitLabels == 1);
n_inh = sum(ctc.UnitLabels == 2);
n_nan = sum(isnan(ctc.UnitLabels));
fprintf('Classification: %i excitatory | %i inhibitory | %i unclassified\n', ...
    n_exc, n_inh, n_nan);

%% 6 — UMAP sanity check  (optional — comment out to skip)

plotUMAPSanityCheck(ctc);

%% 7 — Inspect classification quality

% --- ACG images (replaces the manual imagesc block in clf_all_batches.m) ---
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

% --- I/E ratio per chip ---
[chip_idx, chip_labels] = rg.combineMetadataIndices(ctc.UnitList, "ChipID");
in_count  = histcounts(chip_idx(ctc.UnitLabels == 2));
all_count = histcounts(chip_idx);
ie_per_chip = in_count ./ all_count;
fprintf('Inhibitory fraction per chip: ');
fprintf('%.2f ', ie_per_chip);
fprintf('\n');

%% 8 — E/I network analysis with EIAnalyzer

% old: generate_unit_spk_mat + create_activity_cutout + plot_EI_network_activity
%      + generate_EI_burst_cutouts + calculate_correlations
%      (separate function calls, manual per-culture loop)

eia_params = struct();
eia_params.Activity.BinSize           = 0.01;       % s per bin (10 ms)
eia_params.Activity.SecCutout         = [0, 1200];  % analysis window (s)
eia_params.BurstDetection.Threshold   = 2;          % I/E ratio defining burst
eia_params.BurstDetection.SmoothWindow= 9;          % median filter width (bins)
eia_params.BurstDetection.PeakCutout  = 150;        % bins either side of peak
eia_params.BurstDetection.PeakHeight  = 0.1;
eia_params.BurstDetection.PeakProminence = 0.1;
eia_params.BurstDetection.MinPeakDistance = 6;
eia_params.BurstDetection.SelectedVariant = 4;      % 4 = both clusters aligned

eia = EIAnalyzer(ctc, eia_params);

% Decompose population activity into E and I firing-rate traces + I/E ratio.
% old: create_activity_cutout(binned_mat, bin_size, sec_cutout, responder_ids, non_responder_ids)
eia.computeActivity();

% Detect burst state: bins where I/E < Threshold AND E > noise floor.
% old: z_binary = plot_EI_network_activity(activity, threshold, smooth_window)
eia.detectBursts();

% Align burst cutouts, cluster into 2 types via t-SNE + k-medoids.
% old: burst_mats = generate_EI_burst_cutouts(activity, z_binary, peak_cutout, ph, pp)
eia.extractBurstCutouts();

% Correlate individual unit firing with population I and E traces.
% old: results = calculate_correlations(results)
eia.computeCorrelations();

%% 9 — Visualise

% Plot firing-rate trace coloured by burst state for the first culture.
% old: plot_EI_network_activity(activity, threshold, smooth_window)
eia.PlotNetworkActivity(1);

% Normalised burst shapes across all cultures.
% old: manual per-culture loop + inline norm in clf_all_batches.m lines 196-218
[norm_bursts, norm_ie] = eia.normalizeBurstCutouts();

figure('Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact');
nexttile; imagesc(norm_ie);     title('I/E ratio (normalised)',       'FontWeight', 'normal'); colorbar
nexttile; imagesc(norm_bursts); title('Total activity (normalised)',  'FontWeight', 'normal'); colorbar
nexttile; plot(nanmean(norm_ie)); hold on; plot(nanmean(norm_bursts));
          legend('I/E', 'Total'); xlabel('Bin'); box off
nexttile; plot(gradient(nanmean(norm_bursts)));
          xlabel('Bin'); title('d/dt total activity', 'FontWeight', 'normal'); box off

%% 10 — Save

save_dir  = fullfile(deephys_root, 'Data', 'Classification');
if ~isfolder(save_dir); mkdir(save_dir); end
save_path = fullfile(save_dir, 'cell_type_classification_results');
save(save_path, 'ctc', 'eia');
fprintf('Saved to: %s.mat\n', save_path);
