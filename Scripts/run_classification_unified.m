% run_classification_unified.m
%
% Unified cell-type classification script.  All data-source combinations are
% controlled from §1 (Configuration) — no other section needs editing.
%
% Supported modes  (set train_source / test_source in §1):
%
%   train_source   test_source   Method called
%   ──────────────────────────────────────────────────────────────────────────
%   "DeePhys"    → "DeePhys"    ctc.classifyUnits()
%   "DeePhys"    → "External"   ctc.classifyExternalUnits(wf_ext, acg_ext, sr_ext)
%   "External"   → "DeePhys"    ctc.classifyUnitsWithExternalTrain(wf_ext, acg_ext, y_ext, sr_ext)
%
% Requirements:
%   - DeePhys on path
%   - run_umap toolbox on path (MATLAB File Exchange)
%   - Kilosort outputs (spike_times.npy etc.) for DeePhys recordings

%% 0 — Paths

deephys_root = "/home/phornauer/Git/DeePhysDev/DeePhys";
addpath(genpath(deephys_root));

%% 1 — Configuration  (edit only this section)

% ── Data sources ─────────────────────────────────────────────────────────────
train_source = "DeePhys";   % "DeePhys" | "External"
test_source  = "DeePhys";   % "DeePhys" | "External"

% ── External data ─────────────────────────────────────────────────────────────
% Fill in whichever side(s) are "External".
% ACGs must be computed at Harmonization.ACGBinSize / ACGLag (see below).
% Waveforms at any sampling rate; resampled automatically.
wf_ext  = [];       % (N_samples × N_units) raw unnormalised waveforms
acg_ext = [];       % (N_bins × N_units)    autocorrelograms
sr_ext  = 30000;    % Hz — acquisition sampling rate of external data
y_ext   = [];       % (1 × N_units) labels  1=excitatory  2=inhibitory
                    %   required only when train_source == "External"

% ── DeePhys recording selection ───────────────────────────────────────────────
batch_ids = 1;
week_ids  = 3;
rg_params.Selection.Inclusion = {{'Region', 'Cortex'}, {'AAV', 0, 128, 129}};
rg_params.Selection.Exclusion = {{'Concentration', inf}};

% Metadata field to identify the treatment culture (bootstrap inhibitory candidate
% detection).  Used only when train_source == "DeePhys"; ignored otherwise.
bootstrap_filter = {'AAV', 128};

% ── CellTypeClassifier parameters ─────────────────────────────────────────────
parameters = struct();

% Bootstrap (DeePhys train only)
parameters.Bootstrap.PreCutout  = [0, 1200];    % s — baseline window
parameters.Bootstrap.PostCutout = [6000, 7200]; % s — post-treatment window
parameters.Bootstrap.BinSize    = 20;           % s per bin
parameters.Bootstrap.NIter      = 1000;
parameters.Bootstrap.Alpha      = 1e-10;

% UMAP
parameters.UMAP.UnitFeatures    = ["FullACG", "ReferenceWaveform"];
parameters.UMAP.NormalizationVar= "ChipID";
parameters.UMAP.GroupingVar     = "Concentration";
parameters.UMAP.GroupingValues  = 0;            % baseline recordings only
parameters.UMAP.NNeighbors      = 100;
parameters.UMAP.MinDist         = 0.1;
parameters.UMAP.Spread          = 1;
parameters.UMAP.NDims           = 10;

% Outlier detection (DeePhys train only)
parameters.OutlierDetection.ContaminationFraction = 0.5;
parameters.OutlierDetection.NObsPerLearner        = 50;
parameters.OutlierDetection.DistancePercentile    = 80;

% Harmonization — shared target spec applied to both DeePhys and external data
parameters.Harmonization.ACGBinSize                 = 0.0005;   % s  (0.5 ms)
parameters.Harmonization.ACGLag                     = 0.1;      % s  (100 ms → 401 bins)
parameters.Harmonization.ACGSource                  = "FullACG"; % "FullACG" | "ACG"
parameters.Harmonization.WaveformTargetSamplingRate = 300000;    % Hz

% ── EI analysis parameters (DeePhys test only) ────────────────────────────────
eia_params = struct();
eia_params.Activity.BinSize              = 0.01;
eia_params.Activity.SecCutout            = [0, 1200];
eia_params.BurstDetection.Threshold      = 2;
eia_params.BurstDetection.SmoothWindow   = 9;
eia_params.BurstDetection.PeakCutout     = 150;
eia_params.BurstDetection.PeakHeight     = 0.1;
eia_params.BurstDetection.PeakProminence = 0.1;
eia_params.BurstDetection.MinPeakDistance= 6;
eia_params.BurstDetection.SelectedVariant= 4;

%% 2 — Load DeePhys data

full_rec_array = load_cortical_cultures(batch_ids, week_ids);
rg = RecordingGroup(full_rec_array, rg_params);
rg.Cultures = rg.groupCultures('RecordingDate');

%% 3 — Initialise CellTypeClassifier

ctc = CellTypeClassifier(rg, parameters);

%% 4 — Populate unit list  (always required)
%
% When train_source == "DeePhys": bootstrap test flags inhibitory candidates.
%   bootstrap_filter restricts which cultures are tested (e.g. {'AAV', 128}).
% When train_source == "External": bootstrap results are unused; the empty
%   filter {} collects all units into ctc.UnitList without candidate filtering.

if train_source == "DeePhys"
    ctc.identifyResponsiveUnits(bootstrap_filter);
else
    ctc.identifyResponsiveUnits({});   % populate UnitList only
end

%% 5 — Generate training labels  (DeePhys train only)

if train_source == "DeePhys"
    ctc.generateTrainLabels();
    fprintf('Training set: %i inhibitory + %i excitatory\n', ...
        sum(ctc.TrainLabels.sorted_y_train == 2), ...
        sum(ctc.TrainLabels.sorted_y_train == 1));
end

%% 6 — Classify

if     train_source == "DeePhys"  && test_source == "DeePhys"
    ctc.classifyUnits();
    result_labels = ctc.UnitLabels;

elseif train_source == "DeePhys"  && test_source == "External"
    result_labels = ctc.classifyExternalUnits(wf_ext, acg_ext, sr_ext);

elseif train_source == "External" && test_source == "DeePhys"
    ctc.classifyUnitsWithExternalTrain(wf_ext, acg_ext, y_ext, sr_ext);
    result_labels = ctc.UnitLabels;

else
    error('Unsupported combination: train_source="%s", test_source="%s".', ...
        train_source, test_source);
end

fprintf('Result: %i excitatory | %i inhibitory\n', ...
    sum(result_labels == 1, 'omitnan'), ...
    sum(result_labels == 2, 'omitnan'));

%% 7 — UMAP sanity check

plotUMAPSanityCheck(ctc);

%% 8 — Inspect classification quality  (DeePhys test only)

if test_source == "DeePhys"
    % ACG images — uses whichever ACGSource was selected
    acg_prop  = char(ctc.Parameters.Harmonization.ACGSource);
    all_acgs  = horzcat(ctc.UnitList.(acg_prop))';
    norm_acgs = all_acgs ./ max(all_acgs, [], 2);

    inh_acgs = norm_acgs(ctc.UnitLabels == 2, :);
    exc_acgs = norm_acgs(ctc.UnitLabels == 1, :);

    [sorted_inh, ~] = sortACGsByPeak(inh_acgs, "max");
    [sorted_exc, ~] = sortACGsByPeak(exc_acgs, "max");

    figure('Color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    nexttile; imagesc(sorted_inh); title('Inhibitory ACGs', 'FontWeight', 'normal'); colorbar
    nexttile; imagesc(sorted_exc); title('Excitatory ACGs',  'FontWeight', 'normal'); colorbar

    % Waveforms and ACGs split by cell type; TrainLabels only available for DeePhys train
    train_labels_arg = [];
    if train_source == "DeePhys"
        train_labels_arg = ctc.TrainLabels;
    end
    plotCellTypeFeatures(ctc.UnitList, ctc.UnitLabels, train_labels_arg);

    % I/E fraction per chip
    [chip_idx, ~] = rg.combineMetadataIndices(ctc.UnitList, "ChipID");
    inh_count = histcounts(chip_idx(ctc.UnitLabels == 2));
    all_count = histcounts(chip_idx);
    ie_per_chip = inh_count ./ all_count;
    fprintf('Inhibitory fraction per chip:'); fprintf('  %.2f', ie_per_chip); fprintf('\n');
end

%% 9 — E/I network analysis  (DeePhys test only)

if test_source == "DeePhys"
    eia = EIAnalyzer(ctc, eia_params);
    eia.computeActivity();
    eia.detectBursts();
    eia.extractBurstCutouts();
    eia.computeCorrelations();

    %% 10 — Visualise burst dynamics

    eia.PlotNetworkActivity(1);

    eia.normalizeBurstCutouts();
    norm_bursts = eia.NormalizedCutouts.bursts;
    norm_ie     = eia.NormalizedCutouts.ie;

    figure('Color', 'w');
    tiledlayout(2, 2, 'TileSpacing', 'compact');
    nexttile; imagesc(norm_ie);     title('I/E ratio (normalised)',      'FontWeight', 'normal'); colorbar
    nexttile; imagesc(norm_bursts); title('Total activity (normalised)', 'FontWeight', 'normal'); colorbar
    nexttile; plot(nanmean(norm_ie)); hold on; plot(nanmean(norm_bursts));
              legend('I/E', 'Total'); xlabel('Bin'); box off
    nexttile; plot(gradient(nanmean(norm_bursts)));
              xlabel('Bin'); title('d/dt total activity', 'FontWeight', 'normal'); box off
end

%% 11 — Save

save_dir  = fullfile(deephys_root, 'Data', 'Classification');
if ~isfolder(save_dir); mkdir(save_dir); end

mode_tag  = sprintf('%s_to_%s', lower(train_source), lower(test_source));
save_path = fullfile(save_dir, ['cell_type_classification_' mode_tag]);

if test_source == "DeePhys" && exist('eia', 'var')
    save(save_path, 'ctc', 'eia');
else
    save(save_path, 'ctc');
end
fprintf('Saved to: %s.mat\n', save_path);
