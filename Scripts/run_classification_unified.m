% run_classification_unified.m
%
% Unified classification script — all data-source combinations controlled
% from §1 (Configuration).
%
% Supported modes (set train_source / test_source):
%
%   train_source   test_source   Method
%   ─────────────────────────────────────────────────────────────────────
%   "DeePhys"    → "DeePhys"    ctc.classifyUnits()
%   "DeePhys"    → "External"   ctc.classifyExternalUnits(wf,acg,sr)
%   "External"   → "DeePhys"    ctc.classifyUnitsWithExternalTrain(wf,acg,y,sr)
%
% Pipeline (DeePhys → DeePhys):
%   generateTrainLabels: feature extraction + per-group z-score (cached),
%     global z-score fit on training subset, unsupervised UMAP, outlier detection
%   classifyUnits: reuses feature cache, applies stored NormalizationParams,
%     supervised UMAP, distance-weighted kNN confidence

%% 0 — Setup

deephys_root = 'C:/Users/pjhor/Documents/Code/DeePhys';
addpath(genpath(deephys_root));

%% 1 — Configuration

train_source = "DeePhys";   % "DeePhys" | "External"
test_source  = "DeePhys";   % "DeePhys" | "External"
save_dir     = '/path/to/save';

% DeePhys data
proc_paths = {
    '/path/to/proc1.mat'
    '/path/to/proc2.mat'
};
bootstrap_filter = {'AAV', 128};   % metadata filter for inhibitory candidates

% External data (fill in whichever side(s) are "External")
wf_ext  = [];       % (N_samples x N_units) raw waveforms
acg_ext = [];       % (N_bins x N_units) autocorrelograms
sr_ext  = 30000;    % waveform sampling rate (Hz)
y_ext   = [];       % (1 x N_units) labels: 1=exc, 2=inh  (train_source=="External" only)

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
parameters.UMAP.TargetWeight       = 0.5;
parameters.Harmonization.ACGSource  = 'FullACG';
parameters.Harmonization.ACGBinSize = 0.0005;
parameters.Harmonization.ACGLag     = 0.1;
parameters.OutlierDetection.OutlierAlpha            = 0.01;
parameters.OutlierDetection.DipTestAlpha            = 0.05;
parameters.OutlierDetection.MaxResponsiveComponents = 3;
% CounterexampleRatio: excitatory candidates per inhibitory candidate.
% Default is 1 (balanced). Increase to 2-3 if inhibitory fraction is too high (>30%).
parameters.OutlierDetection.CounterexampleRatio     = 1;

% Optional adaptive parameters (uncomment to enable):
%   parameters.UMAP.AutoNNeighbors              = true;   % n_neighbors = max(15, sqrt(N))
%   parameters.UMAP.AutoConfidenceK             = true;   % kNN k = max(5, sqrt(N_train))
%   parameters.UMAP.FeatureSelection            = true;   % remove low-var / correlated features
%   parameters.Bootstrap.UseFDR                 = true;   % BH correction instead of fixed alpha

%% 2 — Load DeePhys data (when needed)

if train_source == "DeePhys" || test_source == "DeePhys"
    procs = RecordingProcessor.loadMany(proc_paths);
    fs    = FeatureStore.fromProcessors(procs);
    ud    = [];
    for i = 1:numel(procs)
        ud = [ud, procs(i).Units]; %#ok<AGROW>
    end
    ctc = CellTypeClassifier(fs, ud, parameters);
    fprintf('DeePhys: %d units, %d recordings\n', ...
        height(fs.UnitTable), height(fs.RecordingTable));
end

%% 3 — Train on DeePhys (bootstrap + label generation when needed)

if train_source == "DeePhys"
    ctc.identifyResponsiveUnits(bootstrap_filter);
    ctc.generateTrainLabels();
    fprintf('Train set: %d exc + %d inh\n', ...
        sum(ctc.TrainLabels.sorted_y_train == 1), ...
        sum(ctc.TrainLabels.sorted_y_train == 2));
end

%% 4 — Classify

if train_source == "DeePhys" && test_source == "DeePhys"
    ctc.classifyUnits();

elseif train_source == "DeePhys" && test_source == "External"
    assert(~isempty(wf_ext) && ~isempty(acg_ext), ...
        'Fill in wf_ext and acg_ext for External test source.');
    ext_labels = ctc.classifyExternalUnits(wf_ext, acg_ext, sr_ext);
    fprintf('External: %d excitatory, %d inhibitory\n', ...
        sum(ext_labels == 1), sum(ext_labels == 2));

elseif train_source == "External" && test_source == "DeePhys"
    assert(~isempty(wf_ext) && ~isempty(acg_ext) && ~isempty(y_ext), ...
        'Fill in wf_ext, acg_ext, and y_ext for External train source.');
    ctc.classifyUnitsWithExternalTrain(wf_ext, acg_ext, y_ext, sr_ext);
end

if train_source == "DeePhys" || test_source == "DeePhys"
    if ~isempty(ctc.UnitLabels)
        fprintf('Result: %d excitatory | %d inhibitory | %d unclassified\n', ...
            sum(ctc.UnitLabels == 1), sum(ctc.UnitLabels == 2), sum(isnan(ctc.UnitLabels)));
    end
end

%% 5 — E/I analysis (DeePhys test only)

if test_source == "DeePhys" && ~isempty(ctc.UnitLabels)
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
    plotCellTypeFeatures(ctc);
end

%% 6 — Save

if ~isfolder(save_dir); mkdir(save_dir); end
tag       = strjoin([train_source, "to", test_source], "_");
save_path = fullfile(save_dir, "clf_unified_" + tag);
if test_source == "DeePhys" && exist('eia','var')
    save(save_path, 'ctc', 'eia');
else
    save(save_path, 'ctc');
end
fprintf('Saved to %s.mat\n', save_path);
