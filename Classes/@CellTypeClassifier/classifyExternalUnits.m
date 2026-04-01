function ext_labels = classifyExternalUnits(ctc, wf_test, acg_test, sr_test, opts)
% CLASSIFYEXTERNALUNITS  Classify external units using the trained UMAP space.
%
% Builds feature matrices for both DeePhys training units and the external test
% units via buildFeatureMatrix (using ctc.Parameters.Harmonization), then applies
% per-group z-score (DeePhys side only), joint global z-score, and supervised
% UMAP classification.
%
% INPUTS:
%   ctc      - CellTypeClassifier with TrainLabels and UnitDataArray set
%   wf_test  - (N_samples × N_ext) raw waveforms at sr_test
%   acg_test - (N_bins × N_ext)   ACGs; must be computed with the same
%              ACGBinSize / ACGLag as ctc.Parameters.Harmonization
%   sr_test  - waveform sampling rate of the external data in Hz
%
% OPTIONAL NAME-VALUE:
%   ACGBinSize - bin size (s) used to compute acg_test.  When provided,
%                it is checked against ctc.Parameters.Harmonization.ACGBinSize.
%   ACGLag     - half-width (s) used to compute acg_test.  Checked against
%                ctc.Parameters.Harmonization.ACGLag.
%
%   Omitting ACGBinSize / ACGLag issues a warning that provenance is unknown.
%   Use CellTypeClassifier.computeACGs() to generate compatible ACGs from
%   spike times.
%
% Sets:  ctc.Reduction.External  (N_ext × NDims)
% Returns: ext_labels  (1 × N_ext)  1=excitatory, 2=inhibitory

arguments
    ctc      CellTypeClassifier
    wf_test  (:,:) double
    acg_test (:,:) double
    sr_test  (1,1) double
    opts.ACGBinSize double = []
    opts.ACGLag     double = []
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyExternalUnits().');
assert(~isempty(ctc.UnitDataArray), ...
    'ctc.UnitDataArray is empty — classifyExternalUnits requires DeePhys units for training.');

% ── ACG parameter validation ──────────────────────────────────────────────────
ph             = ctc.Parameters.Harmonization;
n_bins_expected = round(ph.ACGLag / ph.ACGBinSize) * 2 + 1;

assert(size(acg_test, 1) == n_bins_expected, ...
    ['acg_test has %d bins but %d are expected ' ...
     '(ACGBinSize=%.4f s, ACGLag=%.4f s).\n' ...
     'Use CellTypeClassifier.computeACGs(spike_times_cell, %.4f, %.4f) to compute compatible ACGs.'], ...
    size(acg_test, 1), n_bins_expected, ph.ACGBinSize, ph.ACGLag, ph.ACGBinSize, ph.ACGLag);

if ~isempty(opts.ACGBinSize) && ~isempty(opts.ACGLag)
    assert(abs(opts.ACGBinSize - ph.ACGBinSize) < 1e-12 && ...
           abs(opts.ACGLag     - ph.ACGLag)     < 1e-12, ...
        ['ACG parameters supplied (BinSize=%.6f, Lag=%.4f) do not match ' ...
         'classifier (BinSize=%.6f, Lag=%.4f).'], ...
        opts.ACGBinSize, opts.ACGLag, ph.ACGBinSize, ph.ACGLag);
elseif isempty(opts.ACGBinSize) || isempty(opts.ACGLag)
    warning('CellTypeClassifier:acgProvenanceUnknown', ...
        ['ACG parameter provenance unknown — verify acg_test was computed with ' ...
         'BinSize=%.4f s, Lag=%.4f s.'], ph.ACGBinSize, ph.ACGLag);
end

labels = ctc.TrainLabels;

% ── Build DeePhys train features ──────────────────────────────────────────────
[wf_dphy, acg_dphy, sr_dphy] = ctc.getOrExtract(ctc.UnitDataArray);
[X_dphy, feat_names]         = buildFeatureMatrix(ctc, wf_dphy, acg_dphy, sr_dphy);

% ── Step 1: per-group normalisation on DeePhys features ──────────────────────
X_dphy(isnan(X_dphy)) = 0;
norm_var = ctc.Parameters.UMAP.NormalizationVar;
if ~isempty(norm_var) && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col  = ctc.FeatureStore.UnitTable.(norm_var);
    [~, ~, iG] = unique(string(group_col), 'stable');
    for g = 1:max(iG)
        mask            = (iG == g);
        X_dphy(mask, :) = normalize(X_dphy(mask, :));
    end
elseif ~isempty(norm_var)
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in UnitTable — per-group z-score skipped.', norm_var);
end
% External test data has no DeePhys metadata so per-group z-score is not
% applied to X_test. Steps 2–4 are handled inside classifyWithFeatureMatrices.

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
