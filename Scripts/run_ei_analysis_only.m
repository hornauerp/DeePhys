% run_ei_analysis_only.m
%
% E/I network analysis from pre-computed classification labels.
%
% Use this script when a CellTypeClassifier has already been run and saved,
% and only the EIAnalyzer step needs to be (re-)run — e.g., to adjust burst
% detection parameters or the analysis window without re-running UMAP.
%
% This recapitulates the inline E/I loop from:
%   U:\Git\Chemogenetics\Scripts\Classification\clf_all_batches.m  (lines 152–224)
%
% Requirements: DeePhys on path, run_umap toolbox on path.

%% 0 — Setup

deephys_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(deephys_root));

%% 1 — Load pre-computed CellTypeClassifier

% Point to the .mat file saved by run_cell_type_classification.m.
save_path = fullfile(deephys_root, 'Data', 'Classification', ...
                     'cell_type_classification_results.mat');

loaded = load(save_path, 'ctc');
ctc    = loaded.ctc;

fprintf('Loaded CellTypeClassifier: %i units (%i inh, %i exc)\n', ...
    length(ctc.UnitLabels), sum(ctc.UnitLabels == 2), sum(ctc.UnitLabels == 1));

%% 2 — Configure EIAnalyzer

% Adjust these parameters without needing to rerun classification.
eia_params = struct();

% Activity decomposition
eia_params.Activity.BinSize        = 0.01;      % s per bin (10 ms)
eia_params.Activity.SecCutout      = [0, 1200]; % analysis window (s)

% Burst detection
eia_params.BurstDetection.Threshold       = 2;    % I/E ratio; lower = stricter
eia_params.BurstDetection.SmoothWindow    = 9;    % median filter width (bins)
eia_params.BurstDetection.PeakCutout      = 150;  % bins either side of burst peak
eia_params.BurstDetection.PeakHeight      = 0.1;
eia_params.BurstDetection.PeakProminence  = 0.1;
eia_params.BurstDetection.MinPeakDistance = 6;    % minimum bins between peaks

% Burst variant:
%   1 = all bursts (cross-corr aligned)
%   2 = cluster 1 only
%   3 = cluster 2 only
%   4 = both clusters mutually aligned  (default, used in paper)
eia_params.BurstDetection.SelectedVariant = 4;

eia = EIAnalyzer(ctc, eia_params);

%% 3 — Run E/I analysis

% Decompose into E/I firing-rate traces
% old: create_activity_cutout(binned_mat, bin_size, sec_cutout, responder_ids, non_responder_ids)
eia.computeActivity();
fprintf('Activity computed for %i cultures\n', length(eia.Activity));

% Detect burst state (binary, median-filtered)
% old: z_binary = plot_EI_network_activity(activity, threshold, smooth_window)
eia.detectBursts();

% Extract burst cutouts and cluster into two types
% old: burst_mats = generate_EI_burst_cutouts(activity, z_binary, peak_cutout, ph, pp)
eia.extractBurstCutouts();

% Unit↔population correlations
% old: results = calculate_correlations(results)
eia.computeCorrelations();

%% 4 — Visualise

% Network activity with burst-state colouring for each culture
n_cultures = length(eia.Activity);
for c = 1:n_cultures
    % old: plot_EI_network_activity(activity, threshold, smooth_window)
    eia.PlotNetworkActivity(c);
    drawnow;
end

% Normalised burst shapes
% old: manual norm loop in clf_all_batches.m lines 196-218
eia.normalizeBurstCutouts();
norm_bursts = eia.NormalizedCutouts.bursts;
norm_ie     = eia.NormalizedCutouts.ie;

figure('Color', 'w');
tiledlayout(2, 3, 'TileSpacing', 'compact');

nexttile; imagesc(norm_ie);
title('I/E ratio — all cultures', 'FontWeight', 'normal'); colorbar; axis tight

nexttile; imagesc(norm_bursts);
title('Total activity — all cultures', 'FontWeight', 'normal'); colorbar; axis tight

nexttile;
plot(nanmean(norm_ie), 'r'); hold on;
plot(nanmean(norm_bursts), 'k');
legend('I/E', 'Total', 'Location', 'northwest'); xlabel('Bin (10 ms)'); box off

% Split into train and test cultures if training_ids is known
% (mirrors clf_transfer.m lines 89–101)
% Uncomment and set train_culture_idx if applicable:
% train_culture_idx = 1:8;
% test_culture_idx  = setdiff(1:n_cultures, train_culture_idx);
% nexttile;
% plot(nanmean(norm_ie(train_culture_idx,:))); hold on
% plot(nanmean(norm_bursts(train_culture_idx,:)));
% title('Training cultures', 'FontWeight', 'normal'); box off
% nexttile;
% plot(nanmean(norm_ie(test_culture_idx,:))); hold on
% plot(nanmean(norm_bursts(test_culture_idx,:)));
% title('Test cultures', 'FontWeight', 'normal'); box off

%% 5 — Save updated EIAnalyzer

save_path_eia = fullfile(deephys_root, 'Data', 'Classification', ...
                         'ei_analysis_results');
save(save_path_eia, 'eia');
fprintf('EIAnalyzer saved to: %s.mat\n', save_path_eia);
