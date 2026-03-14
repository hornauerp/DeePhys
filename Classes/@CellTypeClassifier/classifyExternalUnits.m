function ext_labels = classifyExternalUnits(ctc, wf_test, acg_test, sr_test)
% CLASSIFYEXTERNALUNITS  Classify external units using the trained UMAP space.
%
% Builds feature matrices for both DeePhys training units and the external test
% units via buildFeatureMatrix (using ctc.Parameters.Harmonization), then applies
% joint z-score normalisation and supervised UMAP classification.
%
% To harmonize DeePhys ACGs to a different resolution than the default (0.5 ms /
% 100 ms), change ctc.Parameters.Harmonization.ACGBinSize and .ACGLag before
% calling this function.  Supply acg_test computed with the same parameters.
%
% INPUTS:
%   ctc      - CellTypeClassifier with TrainLabels and UnitList set
%   wf_test  - (N_samples × N_ext) raw waveforms at sr_test
%   acg_test - (N_bins × N_ext)   ACGs; must satisfy
%              N_bins == round(2*Harmonization.ACGLag/Harmonization.ACGBinSize)+1
%   sr_test  - waveform sampling rate of the external data in Hz
%
% Sets:  ctc.Reduction.External  (N_ext × NDims)
% Returns: ext_labels  (1 × N_ext)  1=excitatory, 2=inhibitory

arguments
    ctc     CellTypeClassifier
    wf_test  (:,:) double
    acg_test (:,:) double
    sr_test  (1,1) double
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyExternalUnits().');
assert(~isempty(ctc.UnitList), ...
    'ctc.UnitList is empty — classifyExternalUnits requires DeePhys units for training.');

labels = ctc.TrainLabels;

% ── Build DeePhys train features using harmonization parameters ───────────────
[wf_dphy, acg_dphy, sr_dphy] = extractUnitWaveformsAndACGs(ctc, ctc.UnitList);
[X_dphy, feat_names]         = buildFeatureMatrix(ctc, wf_dphy, acg_dphy, sr_dphy);

% ── Build external test features ─────────────────────────────────────────────
[X_test, feat_names_test] = buildFeatureMatrix(ctc, wf_test, acg_test, sr_test);

assert(isequal(feat_names, feat_names_test), ...
    ['Feature mismatch: DeePhys (%i) vs external (%i) features.\n' ...
     'Ensure external ACGs use the same bin_size/lag as ctc.Parameters.Harmonization.'], ...
    length(feat_names), length(feat_names_test));

% ── Select training subset and classify ──────────────────────────────────────
X_train = X_dphy(labels.umap_train_idx, :);

[ext_labels, ~, ext_reduction] = classifyWithFeatureMatrices(ctc, ...
    X_train, labels.sorted_y_train, X_test, feat_names);

ctc.Reduction.External = ext_reduction;

fprintf('External classification: %i excitatory, %i inhibitory (%i total)\n', ...
    sum(ext_labels == 1), sum(ext_labels == 2), size(X_test, 1));
end
