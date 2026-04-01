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
%   params.Harmonization.ACGLag     = 0.1;        % ±100 ms lag
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
% Runs unsupervised UMAP to embed all units, then uses an isolation forest
% to identify high-density clusters corresponding to excitatory and inhibitory
% populations. Inhibitory candidates (from Step 3) seed the inhibitory cluster.

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

wf  = ctc.HarmonizedWaveforms;   % (N_samples × N_units)
acg = ctc.HarmonizedACGs;         % (N_bins   × N_units)
sr  = ctc.HarmonizedSR;           % target waveform sampling rate

fprintf('Waveform size : %d × %d at %.0f Hz\n', size(wf,1), size(wf,2), sr);
fprintf('ACG size      : %d × %d\n', size(acg,1), size(acg,2));

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
ctc2.Parameters.Harmonization.ACGLag     = 0.05;    % ±50 ms
ctc2.identifyResponsiveUnits();
ctc2.generateTrainLabels();
ctc2.classifyUnits();

%% 10  Metadata-driven labels (skip bootstrap)
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

%% 11  Inspect UMAP embedding

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
