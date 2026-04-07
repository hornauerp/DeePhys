%% Tutorial 6 — Transfer Learning
%
% Demonstrates three transfer scenarios:
%   A. Cross-culture:  train on subset of cultures, apply to the rest
%   B. Cross-dataset:  train on dataset A, classify dataset B
%   C. External train: classify DeePhys units using labelled external data
%
% Waveform/ACG harmonization handles different sampling rates and bin sizes
% so classifiers can be applied across platforms.
%
% Prerequisites:
%   - Tutorial 4 completed on a source dataset
%   - DeePhys on the MATLAB path

%% ============================================================
%% A  Cross-culture transfer
%% ============================================================

%% A1  Load source dataset and split by culture

fs_file   = '/path/to/FeatureStore.mat';
proc_paths = {'/path/to/proc1.mat', '/path/to/proc2.mat'};

fs    = FeatureStore.load(fs_file);
procs = RecordingProcessor.loadMany(proc_paths);

ud = [];
for i = 1:numel(procs)
    ud = [ud, procs(i).Units]; %#ok<AGROW>
end

% Unique cultures in the dataset
culture_ids = FeatureStore.getCultureIDsForUnits( ...
    fs.UnitTable.RecordingID, fs.MetadataTable, ["ChipID","PlatingDate"]);
unique_cultures = unique(culture_ids, 'stable');
fprintf('Total cultures: %d\n', numel(unique_cultures));

%% A2  Subset: use first half as source, second half as target

n_total  = numel(unique_cultures);
n_source = floor(n_total / 2);
source_cultures = unique_cultures(1:n_source);
target_cultures = unique_cultures(n_source+1:end);

% Get recording IDs belonging to source / target cultures
source_rec_ids = fs.MetadataTable.RecordingID( ...
    ismember(FeatureStore.buildCultureIDs(fs.MetadataTable, ["ChipID","PlatingDate"]), source_cultures));
target_rec_ids = fs.MetadataTable.RecordingID( ...
    ismember(FeatureStore.buildCultureIDs(fs.MetadataTable, ["ChipID","PlatingDate"]), target_cultures));

fs_source = fs.subset(source_rec_ids);
fs_target = fs.subset(target_rec_ids);

% Split UnitData accordingly
source_mask = ismember(string({ud.RecordingID})', source_rec_ids);
target_mask = ismember(string({ud.RecordingID})', target_rec_ids);
ud_source = ud(source_mask);
ud_target = ud(target_mask);

%% A3  Train classifier on source cultures

ctc_source = CellTypeClassifier(fs_source, ud_source);
ctc_source.identifyResponsiveUnits();
ctc_source.generateTrainLabels();
ctc_source.classifyUnits();

%% A4  Apply to target cultures
%
% applyTo creates a new CellTypeClassifier for the target dataset,
% normalizes features using the same pipeline, and applies the trained
% UMAP projection. No re-training is performed on the target.

[target_labels, ctc_target] = ctc_source.applyTo(fs_target, ud_target);

n_exc = sum(target_labels == 1);
n_inh = sum(target_labels == 2);
fprintf('Target: %d excitatory, %d inhibitory\n', n_exc, n_inh);

% Run EI analysis on target
eia = EIAnalyzer(ctc_target);
eia.computeActivity();
eia.detectBursts();

%% A5  Transfer with a different normalization variable

[target_labels2, ~] = ctc_source.applyTo(fs_target, ud_target, ...
    'NormalizationVar', 'ChipID');

%% ============================================================
%% B  Cross-dataset transfer
%% ============================================================

%% B1  Load a second dataset (dataset B)

fs_file_B   = '/path/to/FeatureStore_B.mat';
proc_paths_B = {'/path/to/procB1.mat', '/path/to/procB2.mat'};

fs_B    = FeatureStore.load(fs_file_B);
procs_B = RecordingProcessor.loadMany(proc_paths_B);

ud_B = [];
for i = 1:numel(procs_B)
    ud_B = [ud_B, procs_B(i).Units]; %#ok<AGROW>
end

%% B2  Apply source classifier to dataset B

[labels_B, ctc_B] = ctc_source.applyTo(fs_B, ud_B);
fprintf('Dataset B: %d excitatory, %d inhibitory\n', ...
    sum(labels_B == 1), sum(labels_B == 2));

%% ============================================================
%% C  Classify DeePhys units using external labelled training data
%% ============================================================

%% C1  Prepare external data
%
% External data: waveforms and ACGs from a dataset with known ground-truth
% labels (e.g. patch-clamp identified, or optogenetically tagged).
%
% Requirements:
%   wf_ext  : (N_samples × N_ext) raw unnormalized waveforms at sr_ext Hz
%   acg_ext : (N_bins × N_ext)   ACGs computed with the SAME bin_size and lag
%             as ctc.Parameters.Harmonization (default: 0.5 ms / ±100 ms)
%   y_ext   : (1 × N_ext) class labels: 1=excitatory, 2=inhibitory
%   sr_ext  : sampling rate of the external recording in Hz

% Placeholder — replace with actual external data loading
sr_ext  = 30000;
n_ext   = 50;
n_samp  = 60;   % waveform samples in external recording
n_bins  = round(2 * 0.1 / 0.0005) + 1;   % default lag/binsize
wf_ext  = randn(n_samp, n_ext);
acg_ext = rand(n_bins, n_ext);
y_ext   = [ones(1, 25), 2*ones(1, 25)];  % 25 exc + 25 inh

%% C2  Classify DeePhys units with external training data

ctc_ext = CellTypeClassifier(fs, ud);
ctc_ext.classifyUnitsWithExternalTrain(wf_ext, acg_ext, y_ext, sr_ext);

labels_ext = ctc_ext.UnitLabels;
fprintf('External train: %d excitatory, %d inhibitory\n', ...
    sum(labels_ext == 1), sum(labels_ext == 2));

%% C3  Classify external units with DeePhys-trained classifier
%
% Inverse scenario: use DeePhys train labels to classify external test units.

ctc_dphy = CellTypeClassifier(fs, ud);
ctc_dphy.identifyResponsiveUnits();
ctc_dphy.generateTrainLabels();

% External test waveforms/ACGs (different device/sampling rate)
n_test_ext  = 30;
wf_test_ext  = randn(n_samp, n_test_ext);
acg_test_ext = rand(n_bins, n_test_ext);

ext_labels = ctc_dphy.classifyExternalUnits(wf_test_ext, acg_test_ext, sr_ext);
fprintf('Classified external: %d excitatory, %d inhibitory\n', ...
    sum(ext_labels == 1), sum(ext_labels == 2));

%% C4  Waveform harmonization note
%
% If external data has a different sampling rate, buildFeatureMatrix
% (called internally) resamples waveforms to ctc.Parameters.Harmonization.
% WaveformTargetSamplingRate (default 120 kHz) before computing features.
%
% For ACGs, you must supply them pre-computed at the same BinSize and Lag
% as ctc.Parameters.Harmonization.ACGBinSize and .ACGLag.
% Adjust the classifier parameters if needed:
%
%   ctc.Parameters.Harmonization.ACGBinSize = 0.001;  % 1 ms bins
%   ctc.Parameters.Harmonization.ACGLag     = 0.05;   % ±50 ms
%   % then recompute acg_ext at 1 ms bins / ±50 ms before calling classify*
