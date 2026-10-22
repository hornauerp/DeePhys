function ext_labels = classifyExternalUnits(ctc, wf_test, acg_test, sr_test, opts)
% CLASSIFYEXTERNALUNITS  Classify external units via a joint UMAP with DeePhys units.
%
% Builds a joint unsupervised UMAP embedding of all DeePhys units plus the
% external (unlabelled) units. DeePhys training labels from generateTrainLabels
% seed graph label propagation, which then assigns labels to the external units.
%
% Feature space consistency is enforced by applying the stored NormalizationParams
% (from generateTrainLabels) to both DeePhys and external feature matrices before
% embedding, so both populations are placed in the same normalized feature space.
%
% INPUTS:
%   ctc      - CellTypeClassifier with TrainLabels, NormalizationParams, and
%              UnitDataArray set (run identifyResponsiveUnits + generateTrainLabels first)
%   wf_test  - (N_samples × N_ext) raw waveforms at sr_test
%   acg_test - (N_bins × N_ext)   ACGs computed with the same ACGBinSize /
%              ACGLag as ctc.Parameters.Harmonization
%   sr_test  - waveform sampling rate of the external data in Hz
%
% OPTIONAL NAME-VALUE:
%   ACGBinSize        - bin size (s) used to compute acg_test. When provided,
%                       it is checked against ctc.Parameters.Harmonization.ACGBinSize.
%   ACGLag            - half-width (s) used to compute acg_test. Checked against
%                       ctc.Parameters.Harmonization.ACGLag.
%   WaveformPreTrough - ms of waveform captured before the trough. When provided,
%                       it is checked against ctc.Parameters.Harmonization.WaveformPreTrough.
%   WaveformPostTrough - ms of waveform captured after the trough. Checked against
%                        ctc.Parameters.Harmonization.WaveformPostTrough.
%   PerGroupZScore     - "auto" (default) | true | false. DeePhys features are z-scored
%                        per Parameters.UMAP.NormalizationVar group before the global
%                        normalisation; external features must get the same step or they
%                        land on a different scale (an island in the joint graph that label
%                        propagation never reaches). "auto" mirrors the classifier.
%                        Caveat: a per-group z-score centres each group on its own mean, so
%                        it assumes groups of similar class composition. Validation sets
%                        enriched for one class (e.g. optotagged interneurons) are biased by
%                        it; train with NormalizationVar = "" for such transfers.
%   ExternalGroups     - (N_ext x 1) group labels for that z-score (default: one group).
%   Cutoff             - threshold on P(inhibitory): 0.5 (default; the calibrated point within
%                        the training domain) or "auto": inferred from the external scores
%                        alone (two-component mixture on logit P, used only if clearly
%                        bimodal; inferExternalCutoff). Use "auto" when the external data come
%                        from another preparation (e.g. in vivo) and attach to the training
%                        graph through a biased set of units: the ranking then stays good but
%                        all scores shift. The culture labels are not affected. Result in
%                        Reduction.External.Cutoff / CutoffInfo.
%                        A per-group EM prior adjustment (as estimateComposition does for
%                        cultures) is not offered: external scores are shifted, not only
%                        re-weighted, and the EM then drifts towards one class.
%
%   Omitting ACG / waveform provenance parameters issues a warning.
%   Use CellTypeClassifier.computeACGs() to generate compatible ACGs from
%   spike times.
%
% Sets:  ctc.Reduction.External.Embedding    (N_ext × NDims) external embedding
%        ctc.Reduction.External.DeePhysEmbedding (n_unique × NDims) DeePhys embedding
%        ctc.Reduction.External.Confidence   (1 × N_ext) label confidence
%        ctc.Reduction.External.JointGraph   (n_unique+N_ext square) graph used for propagation
%        ctc.Reduction.External.Probability  (N_ext × 2) [P(exc) P(inh)] after propagation
% Returns: ext_labels  (1 × N_ext)  1=excitatory, 2=inhibitory, NaN=not reached by
%          label propagation (no graph path to a DeePhys training unit)

arguments
    ctc      CellTypeClassifier
    wf_test  (:,:) double
    acg_test (:,:) double
    sr_test  (1,1) double
    opts.ACGBinSize         double = []
    opts.ACGLag             double = []
    opts.WaveformPreTrough  double = []
    opts.WaveformPostTrough double = []
    opts.PerGroupZScore            = "auto"
    opts.ExternalGroups            = []
    opts.Cutoff                    = 0.5
end

auto_cutoff = (isstring(opts.Cutoff) || ischar(opts.Cutoff)) && lower(string(opts.Cutoff)) == "auto";
assert(auto_cutoff || (isnumeric(opts.Cutoff) && isscalar(opts.Cutoff) && opts.Cutoff > 0 && opts.Cutoff < 1), ...
    'CellTypeClassifier:badCutoff', 'Cutoff must be a number in (0,1) or "auto".');

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyExternalUnits().');
assert(~isempty(ctc.NormalizationParams), ...
    'NormalizationParams not set — run generateTrainLabels() to store them.');
assert(~isempty(ctc.UnitDataArray), ...
    'ctc.UnitDataArray is empty — classifyExternalUnits requires DeePhys units.');

% ── ACG parameter validation ───────────────────────────────────────────────────
ph              = ctc.Parameters.Harmonization;
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

% ── Waveform window parameter validation ───────────────────────────────────────
if ~isempty(opts.WaveformPreTrough) && ~isempty(opts.WaveformPostTrough)
    assert(abs(opts.WaveformPreTrough  - ph.WaveformPreTrough)  < 1e-9 && ...
           abs(opts.WaveformPostTrough - ph.WaveformPostTrough) < 1e-9, ...
        ['Waveform window supplied (PreTrough=%.3f ms, PostTrough=%.3f ms) does not match ' ...
         'classifier (PreTrough=%.3f ms, PostTrough=%.3f ms). ' ...
         'Waveforms must be captured with at least the training window around the trough.'], ...
        opts.WaveformPreTrough, opts.WaveformPostTrough, ph.WaveformPreTrough, ph.WaveformPostTrough);
else
    warning('CellTypeClassifier:waveformProvenanceUnknown', ...
        ['Waveform window provenance unknown — verify wf_test was captured with at least ' ...
         '%.3f ms before and %.3f ms after the trough.'], ph.WaveformPreTrough, ph.WaveformPostTrough);
end

labels = ctc.TrainLabels;
np     = ctc.NormalizationParams;
p_umap = ctc.Parameters.UMAP;
p_conf = ctc.Parameters.Classification;

% ── DeePhys features: apply stored NormalizationParams ─────────────────────────
% buildNormalizedFeatures is idempotent; result is cached on the instance.
ctc.buildNormalizedFeatures();
nf = ctc.NormalizedFeatures;

X_dphy = normalize(nf.X_pergroup, 'center', np.mu_global, 'scale', np.sigma_global);
X_dphy(:, np.nan_cols) = [];
X_dphy = X_dphy ./ np.scale;
if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
    X_dphy = X_dphy(:, np.feature_selection_mask);
end
feat_names_final = np.feat_names_trimmed;   % feature names after all removals
n_dphy = size(X_dphy, 1);

% ── External features: apply same stored NormalizationParams ───────────────────
% This places external units in the same normalized feature space as DeePhys,
% which is required for the joint UMAP geometry to be meaningful.
if isfield(np, 'acg_sigma_ms') && np.acg_sigma_ms > 0   % same ACG smoothing as the DeePhys units
    acg_test = smoothACG(acg_test, np.acg_sigma_ms / 1000 / ph.ACGBinSize);
end
[X_ext_raw, ~] = buildFeatureMatrix(ctc, wf_test, acg_test, sr_test);
n_ext = size(X_ext_raw, 1);

assert(size(X_ext_raw, 2) == numel(np.mu_global), ...
    ['Feature count mismatch: external test has %d raw features but ' ...
     'NormalizationParams expects %d. Ensure external ACGs use the same ' ...
     'bin size / lag as ctc.Parameters.Harmonization.'], ...
    size(X_ext_raw, 2), numel(np.mu_global));

missing = isnan(X_ext_raw);   % filled with the (group) mean below, as in buildNormalizedFeatures

% Mirror the per-group z-score that buildNormalizedFeatures applied to the DeePhys units.
norm_var   = string(p_umap.NormalizationVar);
uses_group = strlength(norm_var) > 0 && ...
    ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames));
if string(opts.PerGroupZScore) == "auto"
    do_group = uses_group;
else
    do_group = logical(opts.PerGroupZScore);
end
if do_group
    g_ext = opts.ExternalGroups;
    if isempty(g_ext), g_ext = ones(n_ext, 1); end
    assert(numel(g_ext) == n_ext, 'ExternalGroups must have one label per external unit.');
    [~, ~, iG] = unique(string(g_ext(:)), 'stable');
    for g = 1:max(iG)
        X_ext_raw(iG == g, :) = normalize(X_ext_raw(iG == g, :));
    end
    X_ext_raw(isnan(X_ext_raw)) = 0;   % missing values and constant columns within a group
else
    if uses_group
        warning('CellTypeClassifier:externalScaleMismatch', ...
            ['The classifier z-scores DeePhys features per %s group but PerGroupZScore=false: ' ...
             'external features are on a different scale.'], norm_var);
    end
    [~, c] = find(missing);
    X_ext_raw(missing) = np.mu_global(c);   % maps to 0 in the global z-score
end

X_ext = normalize(X_ext_raw, 'center', np.mu_global, 'scale', np.sigma_global);
X_ext(:, np.nan_cols) = [];
X_ext = X_ext ./ np.scale;
if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
    X_ext = X_ext(:, np.feature_selection_mask);
end

% Safety: remove any residual bad columns (can arise when external data has
% feature values outside the range seen during generateTrainLabels).
bad_cols = any(isnan(X_dphy) | isinf(X_dphy), 1) | ...
           any(isnan(X_ext)  | isinf(X_ext),  1);
if any(bad_cols)
    X_dphy(:, bad_cols)         = [];
    X_ext(:, bad_cols)          = [];
    feat_names_final(bad_cols)  = [];
end

% ── Dimension-normalised feature group weighting ───────────────────────────────
X_dphy = applyFeatureGroupWeights(X_dphy, feat_names_final, p_umap);
X_ext  = applyFeatureGroupWeights(X_ext,  feat_names_final, p_umap);

% ── Unsupervised UMAP on combined DeePhys + external set ───────────────────────
X_combined = [X_dphy; X_ext];
N_all = n_dphy + n_ext;

if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_all)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_all);
else
    n_neighbors = p_umap.NNeighbors;
end

fprintf('classifyExternalUnits: joint graph on %d units (%d DeePhys + %d external)\n', ...
    N_all, n_dphy, n_ext);

% ── Joint unit graph (same construction as the main pipeline) ──────────────────
[G, reduction] = jointGraph(ctc, X_combined, n_neighbors, feat_names_final);
if isempty(reduction), reduction = nan(N_all, 0); end
row_sums = sum(G, 2);
row_sums(row_sums == 0) = 1;
W = G ./ row_sums;

% ── Label propagation: DeePhys training units are seeds ────────────────────────
% DeePhys training unit indices (1-based, into unique-unit space) and their labels.
% External units at rows n_dphy+1 : end receive labels via propagation.
train_ids = labels.sorted_train_ids;
y_train   = labels.sorted_y_train;

F_label = ones(N_all, 2) * 0.5;
for ti = 1:numel(train_ids)
    idx = train_ids(ti);
    if y_train(ti) == 1
        F_label(idx, :) = [1, 0];
    else
        F_label(idx, :) = [0, 1];
    end
end

max_iter = p_conf.GraphMaxIter;
tol      = p_conf.GraphConvergenceTol;
for iter = 1:max_iter
    F_new = W * F_label;
    F_row_sum = sum(F_new, 2);
    F_row_sum(F_row_sum == 0) = 1;
    F_new = F_new ./ F_row_sum;
    % Clamp DeePhys training labels
    for ti = 1:numel(train_ids)
        idx = train_ids(ti);
        if y_train(ti) == 1
            F_new(idx, :) = [1, 0];
        else
            F_new(idx, :) = [0, 1];
        end
    end
    if max(abs(F_new(:) - F_label(:))) < tol
        fprintf('Graph label propagation converged at iteration %d\n', iter);
        F_label = F_new;
        break
    end
    F_label = F_new;
end
if iter == max_iter
    fprintf('Graph label propagation reached max iterations (%d)\n', max_iter);
end

% Units without a graph path to any seed keep the uninformed [0.5 0.5]; flag them
% before normalisation (previously the tie rule silently called them excitatory).
unreached = all(abs(F_label - 0.5) < 1e-12, 2);

% Seed-count normalisation, as in classifyUnits (no-op for balanced training sets).
n_seed = [sum(y_train == 1), sum(y_train == 2)];
if all(n_seed > 0) && n_seed(1) ~= n_seed(2)
    F_label = F_label ./ n_seed;
    F_label = F_label ./ max(sum(F_label, 2), eps);
end

% ── Extract external unit labels and confidence ────────────────────────────────
ext_labels     = nan(1, n_ext);
ext_confidence = nan(1, n_ext);
n_unreached    = sum(unreached(n_dphy+1:end));
if n_unreached > 0
    warning('CellTypeClassifier:externalUnreached', ...
        ['%d/%d external units have no graph path to a DeePhys training unit and are left ' ...
         'unlabelled (NaN). Check that external features are normalised like the DeePhys ones.'], ...
        n_unreached, n_ext);
end
% Decision threshold on P(inh): 0.5 (calibrated within the training domain) or inferred
% from the external scores alone (Cutoff = "auto", see inferExternalCutoff).
p_inh_ext = F_label(n_dphy+1:end, 2);
p_inh_ext(unreached(n_dphy+1:end)) = NaN;
cut_info = struct('method', "fixed", 'dBIC', NaN, 'inhibitory_share', NaN, 'bimodal', false);
if auto_cutoff
    [cutoff, cut_info] = inferExternalCutoff(p_inh_ext);
    if ~cut_info.bimodal
        warning('CellTypeClassifier:cutoffNotInferred', ...
            'External P(inh) is not bimodal (dBIC %.1f); keeping the cut-off at 0.5.', cut_info.dBIC);
    end
    fprintf(['External cut-off: %.3f (%s; dBIC %.0f, estimated inhibitory share %.2f). ' ...
             'Ranking is unchanged; only the operating point moves.\n'], ...
        cutoff, cut_info.method, cut_info.dBIC, cut_info.inhibitory_share);
else
    cutoff = opts.Cutoff;
end
for i = 1:n_ext
    gi = n_dphy + i;
    if unreached(gi), continue; end
    p_exc = F_label(gi, 1);
    p_inh = F_label(gi, 2);
    if p_inh > cutoff
        ext_labels(i)     = 2;
        ext_confidence(i) = p_inh;
    else
        ext_labels(i)     = 1;
        ext_confidence(i) = p_exc;
    end
end

ctc.Reduction.External = struct( ...
    'Embedding',          reduction(n_dphy+1:end, :), ...
    'DeePhysEmbedding',   reduction(1:n_dphy, :), ...
    'Confidence',         ext_confidence, ...
    'JointGraph',         G, ...          % rows 1:NumDeePhys = DeePhys units, then external units
    'NumDeePhys',         n_dphy, ...
    'Probability',        F_label(n_dphy+1:end, :), ...   % [P(exc) P(inh)] per external unit
    'Cutoff',             cutoff, ...
    'CutoffInfo',         cut_info);

fprintf('External classification: %d excitatory, %d inhibitory, %d unreached (%d total)\n', ...
    sum(ext_labels == 1), sum(ext_labels == 2), n_unreached, n_ext);
end
