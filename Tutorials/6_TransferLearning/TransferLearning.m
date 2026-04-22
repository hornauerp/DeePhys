%% Tutorial 6 — External Data Interoperability
%
% Covers two scenarios where labelled or unlabelled data from a different
% platform or experiment must be classified jointly with DeePhys units:
%
%   C2. External labels (patch-clamp / optogenetics) train a classifier
%       for DeePhys units — no drug-response bootstrap needed.
%   C3. DeePhys-derived training labels classify units from an external
%       recording (different device, sampling rate, etc.).
%
% Both methods build a joint unsupervised UMAP embedding of DeePhys units
% and external units, then apply graph label propagation — the same mechanism
% used by the standard classifyUnits() pipeline.
%
% For cross-culture or cross-dataset classification within DeePhys, use the
% standard pipeline (Tutorial 4): include all units in one CellTypeClassifier
% — generateTrainLabels() and classifyUnits() operate transductively on the
% full NxN UMAP graph without any separate transfer step.
%
% Prerequisites:
%   - Tutorial 1 completed (FeatureStore and RecordingProcessors saved)
%   - DeePhys on the MATLAB path

%% ============================================================
%% C  Classify DeePhys units using external labelled training data
%% ============================================================
%
% Both methods below build a JOINT unsupervised UMAP that embeds DeePhys units
% and external units together in a single shared space, then use graph label
% propagation — the same mechanism as classifyUnits() — to assign labels.
%
% This is stricter than a projection-based approach: both populations must be
% comparable in normalized feature space, so matching ACG parameters and
% waveform harmonization are prerequisites.

%% C1  Prepare external data
%
% External data: waveforms and ACGs from a dataset with known ground-truth
% labels (e.g. patch-clamp identified, or optogenetically tagged).
%
% Requirements:
%   wf_ext  : (N_samples × N_ext) raw unnormalized waveforms at sr_ext Hz
%   acg_ext : (N_bins × N_ext)   ACGs computed with the SAME ACGBinSize and
%             ACGLag as ctc.Parameters.Harmonization (default: 0.5 ms / ±100 ms)
%   y_ext   : (1 × N_ext) class labels: 1=excitatory, 2=inhibitory
%   sr_ext  : sampling rate of the external recording in Hz
%
% Number of ACG bins must equal:  round(ACGLag / ACGBinSize) * 2 + 1
% Use CellTypeClassifier.computeACGs(spike_times, ACGBinSize, ACGLag) if you
% need to compute ACGs from spike times with the correct parameters.

% Placeholder — replace with actual external data loading
sr_ext  = 30000;
n_ext   = 50;
n_samp  = 60;   % waveform samples in external recording
n_bins  = round(2 * 0.1 / 0.0005) + 1;   % default lag/binsize → 401 bins
wf_ext  = randn(n_samp, n_ext);
acg_ext = rand(n_bins, n_ext);
y_ext   = [ones(1, 25), 2*ones(1, 25)];  % 25 exc + 25 inh

%% C2  Classify DeePhys units with external training data
%
% Pipeline:
%   1. buildNormalizedFeatures — per-group z-score on DeePhys units.
%   2. Stack DeePhys (per-group normalized) + external (no per-group normalization).
%   3. Joint global z-score + max(abs) scaling + feature group weighting.
%   4. Unsupervised UMAP on the combined set.
%   5. Graph label propagation: external labels are seeds, DeePhys units receive labels.
%
% Note: identifyResponsiveUnits / generateTrainLabels are NOT required for this path.

ctc_ext = CellTypeClassifier(fs, ud);
ctc_ext.classifyUnitsWithExternalTrain(wf_ext, acg_ext, y_ext, sr_ext);

labels_ext = ctc_ext.UnitLabels;
fprintf('External train: %d excitatory, %d inhibitory\n', ...
    sum(labels_ext == 1, 'omitnan'), sum(labels_ext == 2, 'omitnan'));

% Inspect the joint embedding — DeePhys and external units share the same space
emb_dphy = ctc_ext.Reduction.ExternalTrain.DeePhys;    % (n_unique × NDims)
emb_ext  = ctc_ext.Reduction.ExternalTrain.External;   % (N_train  × NDims)
figure;
scatter(emb_dphy(:,1), emb_dphy(:,2), 10, ctc_ext.UnitLabels( ...
    ctc_ext.NormalizedFeatures.unique_to_rep), 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(emb_ext(:,1), emb_ext(:,2), 30, y_ext', 'filled', ...
    'Marker', '^', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
title('Joint UMAP: DeePhys (circles) + external seeds (triangles)');

%% C3  Classify external units with DeePhys-trained classifier
%
% Inverse scenario: use DeePhys train labels to classify external test units.
%
% Pipeline:
%   1. Apply stored NormalizationParams to both DeePhys and external features
%      (ensures both sides are in the same normalized feature space).
%   2. Unsupervised UMAP on the combined set.
%   3. Graph label propagation: DeePhys training labels are seeds, external
%      units receive labels.
%
% Requires generateTrainLabels() to have been run so NormalizationParams are set.

ctc_dphy = CellTypeClassifier(fs, ud);
ctc_dphy.identifyResponsiveUnits();
ctc_dphy.generateTrainLabels();

% External test waveforms/ACGs (different device/sampling rate)
n_test_ext   = 30;
wf_test_ext  = randn(n_samp, n_test_ext);
acg_test_ext = rand(n_bins, n_test_ext);

ext_labels = ctc_dphy.classifyExternalUnits(wf_test_ext, acg_test_ext, sr_ext);
fprintf('Classified external: %d excitatory, %d inhibitory\n', ...
    sum(ext_labels == 1), sum(ext_labels == 2));

% Confidence and embedding for the external units
ext_conf = ctc_dphy.Reduction.External.Confidence;
emb_ext2 = ctc_dphy.Reduction.External.Embedding;   % (N_ext × NDims)

%% C4  Waveform harmonization note
%
% If external data has a different sampling rate, buildFeatureMatrix
% (called internally) resamples waveforms to
% ctc.Parameters.Harmonization.WaveformTargetSamplingRate (default 120 kHz)
% before computing features — no preprocessing needed on your side.
%
% For ACGs, you must supply them pre-computed at the same BinSize and Lag
% as ctc.Parameters.Harmonization.ACGBinSize and .ACGLag.
% Adjust the classifier parameters if needed:
%
%   ctc.Parameters.Harmonization.ACGBinSize = 0.001;  % 1 ms bins
%   ctc.Parameters.Harmonization.ACGLag     = 0.05;   % ±50 ms
%   % then recompute acg_ext at 1 ms bins / ±50 ms before calling classify*
%
% A mismatch in ACG dimensions triggers an assertion with the expected bin
% count and a suggestion to use CellTypeClassifier.computeACGs().
