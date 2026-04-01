function classifyUnitsWithExternalTrain(ctc, wf_train, acg_train, y_train, sr_train, opts)
% CLASSIFYUNITSWITHEXTERNALTRAIN  Classify DeePhys units trained on external data.
%
% Builds a feature matrix from the supplied external training data, builds a
% second matrix from ctc.UnitDataArray (via the harmonized DeePhys path), applies
% per-group z-score to the DeePhys side, then runs joint global z-score and
% supervised UMAP classification.
%
% Use this when labelled ground-truth data from another source (e.g. a patched
% or optogenetically identified dataset) should serve as the training reference
% for classifying your DeePhys units.
%
% INPUTS:
%   ctc       - CellTypeClassifier with UnitDataArray set
%   wf_train  - (N_samples × N_train) raw unnormalised waveforms at sr_train
%   acg_train - (N_bins × N_train)   ACGs; must be computed with the same
%               ACGBinSize / ACGLag as ctc.Parameters.Harmonization
%   y_train   - (1 × N_train) class labels (1=excitatory, 2=inhibitory)
%   sr_train  - waveform sampling rate of the external training data in Hz
%
% OPTIONAL NAME-VALUE:
%   ACGBinSize - bin size (s) used to compute acg_train.  When provided,
%                it is checked against ctc.Parameters.Harmonization.ACGBinSize.
%   ACGLag     - half-width (s) used to compute acg_train.  Checked against
%                ctc.Parameters.Harmonization.ACGLag.
%
%   Omitting ACGBinSize / ACGLag issues a warning that provenance is unknown.
%   Use CellTypeClassifier.computeACGs() to generate compatible ACGs from
%   spike times.
%
% Sets:  ctc.UnitLabels  (1 = excitatory, 2 = inhibitory)
%        ctc.Reduction.Train  (N_train × NDims)
%        ctc.Reduction.Test   (N_units × NDims)

arguments
    ctc       CellTypeClassifier
    wf_train  (:,:) double
    acg_train (:,:) double
    y_train   (1,:) double
    sr_train  (1,1) double
    opts.ACGBinSize double = []
    opts.ACGLag     double = []
end

assert(~isempty(ctc.UnitDataArray), ...
    'ctc.UnitDataArray is empty — set UnitDataArray before calling classifyUnitsWithExternalTrain().');

% ── ACG parameter validation ──────────────────────────────────────────────────
ph              = ctc.Parameters.Harmonization;
n_bins_expected = round(ph.ACGLag / ph.ACGBinSize) * 2 + 1;

assert(size(acg_train, 1) == n_bins_expected, ...
    ['acg_train has %d bins but %d are expected ' ...
     '(ACGBinSize=%.4f s, ACGLag=%.4f s).\n' ...
     'Use CellTypeClassifier.computeACGs(spike_times_cell, %.4f, %.4f) to compute compatible ACGs.'], ...
    size(acg_train, 1), n_bins_expected, ph.ACGBinSize, ph.ACGLag, ph.ACGBinSize, ph.ACGLag);

if ~isempty(opts.ACGBinSize) && ~isempty(opts.ACGLag)
    assert(abs(opts.ACGBinSize - ph.ACGBinSize) < 1e-12 && ...
           abs(opts.ACGLag     - ph.ACGLag)     < 1e-12, ...
        ['ACG parameters supplied (BinSize=%.6f, Lag=%.4f) do not match ' ...
         'classifier (BinSize=%.6f, Lag=%.4f).'], ...
        opts.ACGBinSize, opts.ACGLag, ph.ACGBinSize, ph.ACGLag);
elseif isempty(opts.ACGBinSize) || isempty(opts.ACGLag)
    warning('CellTypeClassifier:acgProvenanceUnknown', ...
        ['ACG parameter provenance unknown — verify acg_train was computed with ' ...
         'BinSize=%.4f s, Lag=%.4f s.'], ph.ACGBinSize, ph.ACGLag);
end

% ── Build external train features ────────────────────────────────────────────
[X_train, feat_names_train] = buildFeatureMatrix(ctc, wf_train, acg_train, sr_train);

assert(size(X_train, 1) == numel(y_train), ...
    'y_train length (%i) does not match number of external training units (%i).', ...
    numel(y_train), size(X_train, 1));

% ── Build DeePhys test features ───────────────────────────────────────────────
[wf_dphy, acg_dphy, sr_dphy]                      = ctc.getOrExtract(ctc.UnitDataArray);
[X_test, feat_names_test, aligned_wf, norm_acgs]  = buildFeatureMatrix(ctc, wf_dphy, acg_dphy, sr_dphy);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf;
ctc.HarmonizedACGs      = norm_acgs;
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% ── Step 1: per-group normalisation on DeePhys test features ─────────────────
X_test(isnan(X_test)) = 0;
norm_var = ctc.Parameters.UMAP.NormalizationVar;
if ~isempty(norm_var) && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col  = ctc.FeatureStore.UnitTable.(norm_var);
    [~, ~, iG] = unique(string(group_col), 'stable');
    for g = 1:max(iG)
        mask            = (iG == g);
        X_test(mask, :) = normalize(X_test(mask, :));
    end
elseif ~isempty(norm_var)
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in UnitTable — per-group z-score skipped.', norm_var);
end
% External training data has no DeePhys metadata so per-group z-score is not
% applied to X_train. Steps 2–4 are handled inside classifyWithFeatureMatrices.

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
