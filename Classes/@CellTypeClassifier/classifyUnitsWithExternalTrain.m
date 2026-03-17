function classifyUnitsWithExternalTrain(ctc, wf_train, acg_train, y_train, sr_train)
% CLASSIFYUNITSWITHEXTERNALTRAIN  Classify DeePhys units trained on external data.
%
% Builds a feature matrix from the supplied external training data, builds a
% second matrix from ctc.UnitList (via the harmonized DeePhys path), then runs
% joint z-score normalisation and supervised UMAP classification.
%
% Use this when labelled ground-truth data from another source (e.g. a patched
% or optogenetically identified dataset) should serve as the training reference
% for classifying your DeePhys units.
%
% INPUTS:
%   ctc      - CellTypeClassifier with UnitList set
%   wf_train - (N_samples × N_train) raw unnormalised waveforms at sr_train
%   acg_train - (N_bins × N_train)   ACGs; must satisfy
%               N_bins == round(2*Harmonization.ACGLag/Harmonization.ACGBinSize)+1
%   y_train  - (1 × N_train) class labels (1=excitatory, 2=inhibitory)
%   sr_train - waveform sampling rate of the external training data in Hz
%
% Sets:  ctc.UnitLabels  (1 = excitatory, 2 = inhibitory)
%        ctc.Reduction.Train  (N_train × NDims)
%        ctc.Reduction.Test   (N_units × NDims)

arguments
    ctc      CellTypeClassifier
    wf_train  (:,:) double
    acg_train (:,:) double
    y_train   (1,:) double
    sr_train  (1,1) double
end

assert(~isempty(ctc.UnitList), ...
    'ctc.UnitList is empty — set UnitList before calling classifyUnitsWithExternalTrain().');

% ── Build external train features ────────────────────────────────────────────
[X_train, feat_names_train] = buildFeatureMatrix(ctc, wf_train, acg_train, sr_train);

assert(size(X_train, 1) == numel(y_train), ...
    'y_train length (%i) does not match number of external training units (%i).', ...
    numel(y_train), size(X_train, 1));

% ── Build DeePhys test features ───────────────────────────────────────────────
[wf_dphy, acg_dphy, sr_dphy]                      = ctc.getOrExtract(ctc.UnitList);
[X_test, feat_names_test, aligned_wf, norm_acgs]  = buildFeatureMatrix(ctc, wf_dphy, acg_dphy, sr_dphy);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf;
ctc.HarmonizedACGs      = norm_acgs;
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

assert(isequal(feat_names_train, feat_names_test), ...
    ['Feature mismatch: external train (%i) vs DeePhys test (%i) features.\n' ...
     'Ensure external ACGs use the same bin size/lag as ctc.Parameters.Harmonization.'], ...
    length(feat_names_train), length(feat_names_test));

% ── Classify ─────────────────────────────────────────────────────────────────
[Y_pred, train_reduction, test_reduction] = classifyWithFeatureMatrices(ctc, ...
    X_train, y_train, X_test, feat_names_train);

ctc.Reduction.Train = train_reduction;
ctc.Reduction.Test  = test_reduction;
ctc.UnitLabels      = Y_pred;

fprintf('Classified %i excitatory, %i inhibitory DeePhys units (external train, n=%i)\n', ...
    sum(Y_pred == 1), sum(Y_pred == 2), size(X_train, 1));
end
