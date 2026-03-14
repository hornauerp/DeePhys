% run_clf_transfer.m
%
% Cell-type classification and E/I network analysis — transfer-learning layout.
%
% Mirrors the structure of U:\Git\Chemogenetics\Scripts\Classification\clf_transfer.m
% using the new @CellTypeClassifier and @EIAnalyzer infrastructure.
%
% Key logic carried over from clf_transfer.m:
%   1.  Checkpoint guard — section 1 is skipped on re-runs if ctc already exists.
%   2.  Train/test culture split — training_culture_idx selects which recording-date
%       batch was the classifier's training reference (mirrors umap_params.training_ids).
%   3.  5-panel burst visualisation — heatmaps and mean traces for all, train, and
%       test cultures separately (clf_transfer.m lines 89–101).
%
% Parameters are taken from run_cell_type_classification.m (established for this dataset).
%
% Requirements:
%   - DeePhys on path (addpath lines below)
%   - run_umap toolbox on path
%   - Chemogenetics/Helpers on path (load_cortical_cultures)

%% 0 — Setup

deephys_root = "/home/phornauer/Git/DeePhysDev/DeePhys";
addpath(genpath(deephys_root));
addpath("/home/phornauer/Git/Chemogenetics/Helpers");

%% 1 — Load recordings and classify  (guarded — skipped on re-runs)
%
% old: if ~exist("interneuron_idx","var") ... end  (clf_transfer.m lines 1–15)

if ~exist("ctc", "var")

    % --- Load ---
    % old: load_cortical_cultures(batch_ids, week_ids)
    batch_ids = [1, 1, 2, 2];
    week_ids  = [2, 3, 1, 2];
    full_rec_array = load_cortical_cultures(batch_ids, week_ids);

    % Filter recordings: cortex, AAV 0 / 128 / 129, drop wash-out concentration.
    % old: rg_params setup in clf_all_batches.m lines 11–14
    rg_params.Selection.Inclusion = {{'Region', 'Cortex'}, {'AAV', 0, 128, 129}};
    rg_params.Selection.Exclusion = {{'Concentration', inf}};
    rg = RecordingGroup(full_rec_array, rg_params);

    % Group recordings into cultures by recording date.
    % old: dr_rg.Cultures = dr_rg.groupCultures('RecordingDate')
    % new: rg.Cultures is a Culture array (rg.Cultures(c), not {c})
    rg.Cultures = rg.groupCultures('RecordingDate');

    % --- CellTypeClassifier parameters ---
    parameters = struct();

    % Bootstrap firing-rate test (identifies inhibitory candidates)
    % old: pre_cutout, post_cutout, fr_bin_size, n_iter, alpha in clf_all_batches.m lines 27–31
    parameters.Bootstrap.PreCutout  = [0, 1200];    % baseline window (s)
    parameters.Bootstrap.PostCutout = [6000, 7200]; % post-DCZ window (s)
    parameters.Bootstrap.BinSize    = 20;           % bin width (s)
    parameters.Bootstrap.NIter      = 1000;
    parameters.Bootstrap.Alpha      = 1e-10;

    % UMAP embedding
    % old: umap_params in clf_transfer.m lines 18–22
    parameters.UMAP.UnitFeatures     = ["FullACG", "ReferenceWaveform"];
    parameters.UMAP.NormalizationVar = "ChipID";
    parameters.UMAP.GroupingVar      = "Concentration";
    parameters.UMAP.GroupingValues   = 0;           % baseline recordings only
    parameters.UMAP.NNeighbors       = 100;
    parameters.UMAP.MinDist          = 0.1;
    parameters.UMAP.Spread           = 1;
    parameters.UMAP.NDims            = 10;

    % Isolation-forest outlier removal
    % old: outlier_params in clf_transfer.m lines 25–27
    parameters.OutlierDetection.ContaminationFraction = 0.5;
    parameters.OutlierDetection.NObsPerLearner        = 50;
    parameters.OutlierDetection.DistancePercentile    = 80;

    ctc = CellTypeClassifier(rg, parameters);

    % --- Identify inhibitory candidates ---
    % old: [interneuron_idx, unit_list] = identify_interneurons(dr_rg)
    ctc.identifyResponsiveUnits({'AAV', 128});
    fprintf('Responsive (inhibitory candidate) units: %i\n', sum(ctc.ResponsiveUnitIdx));

    % --- Generate training labels ---
    % old: [labels, reduction] = generate_train_labels(dr_rg, unit_list, interneuron_idx, umap_params, outlier_params)
    ctc.generateTrainLabels();

    % --- Classify all units ---
    % old: full_labels = classify_units(dr_rg, unit_list, labels, umap_params)
    ctc.classifyUnits();

    n_exc = sum(ctc.UnitLabels == 1);
    n_inh = sum(ctc.UnitLabels == 2);
    n_nan = sum(isnan(ctc.UnitLabels));
    fprintf('Classification: %i excitatory | %i inhibitory | %i unclassified\n', ...
        n_exc, n_inh, n_nan);

end  % end checkpoint guard

%% 2 — UMAP sanity check  (optional — comment out to skip)

plotUMAPSanityCheck(ctc);

%% 3 — Train / test culture split
%
% old: umap_params.training_ids = 4
%      [train_ids, test_ids] = sort_train_cultures(dr_rg, unit_list, umap_params.training_ids)
%
% training_culture_idx: index/indices into rg.Cultures (recording-date groups) whose
% units were used to train the supervised classifier.
% Adjust to match the batch used during generateTrainLabels.

training_culture_idx = 4;
n_cultures           = length(rg.Cultures);
test_culture_idx     = setdiff(1:n_cultures, training_culture_idx);

fprintf('Training cultures: %s  |  Test cultures: %s\n', ...
    mat2str(training_culture_idx), mat2str(test_culture_idx));

%% 4 — I/E distribution across chips
%
% old: ratio = histcounts(cluster_idx(full_labels==2)) ./ histcounts(cluster_idx)

[chip_idx, ~] = rg.combineMetadataIndices(ctc.UnitList, "ChipID");
in_count  = histcounts(chip_idx(ctc.UnitLabels == 2));
all_count = histcounts(chip_idx);
ie_per_chip = in_count ./ all_count;
fprintf('Inhibitory fraction per chip: ');
fprintf('%.2f ', ie_per_chip);
fprintf('\n');

%% 5 — ACG inspection
%
% old: clf_all_batches.m lines 128–150

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

%% 6 — Validation against ground truth  (optional — comment out to skip)
%
% Compare predicted test-set labels against an independent ground truth.
% Units with NaN ground truth are silently excluded.

gt_field     = 'CellType';     % <-- metadata field carrying ground-truth labels
gt_exc_value = 'PYR';          % <-- value indicating excitatory (class 1)
gt_inh_value = 'PV';           % <-- value indicating inhibitory (class 2)

gt_meta   = arrayfun(@(u) u.MEArecording.Metadata.(gt_field), ...
    ctc.UnitList, 'UniformOutput', false);
gt_labels = nan(1, length(ctc.UnitList));
gt_labels(cellfun(@(v) isequal(v, gt_exc_value), gt_meta)) = 1;
gt_labels(cellfun(@(v) isequal(v, gt_inh_value), gt_meta)) = 2;

% Evaluate on test units only (training labels are known by construction)
test_mask = ctc.TrainLabels.umap_test_idx;
metrics   = validateClassification(ctc.UnitLabels(test_mask), gt_labels(test_mask));

%% 7 — E/I network analysis
%
% old: per-culture loop in clf_all_batches.m / clf_b1_w3.m lines 160–197
%      (generate_unit_spk_mat → create_activity_cutout → plot_EI_network_activity
%       → generate_EI_burst_cutouts → calculate_correlations)

eia_params = struct();
eia_params.Activity.BinSize               = 0.01;       % s per bin (10 ms)
eia_params.Activity.SecCutout             = [0, 1200];  % analysis window (s)
eia_params.BurstDetection.Threshold       = 2;          % I/E ratio threshold
eia_params.BurstDetection.SmoothWindow    = 9;          % median filter width (bins)
eia_params.BurstDetection.PeakCutout      = 150;        % bins either side of burst peak
eia_params.BurstDetection.PeakHeight      = 0.1;
eia_params.BurstDetection.PeakProminence  = 0.1;
eia_params.BurstDetection.MinPeakDistance = 6;
eia_params.BurstDetection.SelectedVariant = 4;          % both clusters aligned

eia = EIAnalyzer(ctc, eia_params);

% old: binned_mat = generate_unit_spk_mat(...)
%      activity   = create_activity_cutout(...)
eia.computeActivity();

% old: z_binary = plot_EI_network_activity(activity, threshold, smooth_window)
eia.detectBursts();

% old: burst_mats = generate_EI_burst_cutouts(activity, z_binary, peak_cutout, ph, pp)
eia.extractBurstCutouts();

% old: results = calculate_correlations(results)
eia.computeCorrelations();

%% 8 — Normalise burst cutouts
%
% old: clf_transfer.m lines 74–87 (mean per culture, row-normalise, extract burst_cutout)
%      new: normalizeBurstCutouts stores the peak-to-end half-window, row-normalised, in eia.NormalizedCutouts

eia.normalizeBurstCutouts();
norm_bursts = eia.NormalizedCutouts.bursts;
norm_ie     = eia.NormalizedCutouts.ie;

%% 9 — Visualise (train / test split)
%
% old: clf_transfer.m lines 89–101
%   imagesc(norm_ie([train_ids;test_ids],:))
%   plot(nanmean(norm_ie(train_ids,:))) / plot(nanmean(norm_ie(test_ids,:)))

sorted_idx = [training_culture_idx(:)', test_culture_idx(:)'];

figure('Color', 'w');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(norm_ie(sorted_idx, :)); colorbar; axis tight
title('I/E ratio  (train | test)', 'FontWeight', 'normal')

nexttile;
imagesc(norm_bursts(sorted_idx, :)); colorbar; axis tight
title('Total activity  (train | test)', 'FontWeight', 'normal')

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

%% 10 — Classify external units  (optional — fill in wf_ext / acg_ext to enable)
%
% Classify units from outside the DeePhys pipeline (e.g. sorted by another tool
% or acquired at a different sampling rate). ACGs must use the same bin size
% and lag as ctc.Parameters.Harmonization; waveforms are resampled automatically.
%
% ACG requirements:
%   bin_size = ctc.Parameters.Harmonization.ACGBinSize  (default: 0.5 ms)
%   lag      = ctc.Parameters.Harmonization.ACGLag       (default: 100 ms → 401 bins)
%
% Waveform requirements:
%   Any sampling rate — set sr_ext to your acquisition rate; resampled automatically.
%   Window should cover the same duration as DeePhys reference waveforms (~2.73 ms).

% wf_ext  = [];    % (N_samples × N_units) unnormalised waveforms
% acg_ext = [];    % (N_bins × N_units)   autocorrelograms (401 bins)
% sr_ext  = 30000; % Hz — set to your acquisition sampling rate
%
% ext_labels = ctc.classifyExternalUnits(wf_ext, acg_ext, sr_ext);
%
% figure('Color','w');
% tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
% nexttile; scatter(ctc.Reduction.External(:,1), ctc.Reduction.External(:,2), ...
%     6, ext_labels, 'filled'); colorbar;
%     xlabel('UMAP 1'); ylabel('UMAP 2'); title('External units — UMAP', 'FontWeight','normal');
% nexttile; histogram(categorical(ext_labels,[1 2],{'Excitatory','Inhibitory'}));
%     ylabel('Count'); title('External unit types', 'FontWeight','normal'); box off

%% 11 — Save
%
% old: save_results = rmfield(results,"binned_mat")  % strip large field
%      save(save_path, "save_results", ...)
% new: eia no longer stores BinnedMats (fixed), safe to save directly

save_dir = fullfile(deephys_root, 'Data', 'Classification');
if ~isfolder(save_dir); mkdir(save_dir); end
save_path = fullfile(save_dir, ...
    strjoin(["clf_transfer", string(training_culture_idx)], "_"));
save(save_path, 'ctc', 'eia', 'norm_bursts', 'norm_ie', ...
     'training_culture_idx', 'test_culture_idx');
fprintf('Saved to %s.mat\n', save_path);
