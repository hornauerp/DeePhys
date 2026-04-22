function classifyUnitsWithExternalTrain(ctc, wf_train, acg_train, y_train, sr_train, opts)
% CLASSIFYUNITSWITHEXTERNALTRAIN  Classify DeePhys units trained on external data.
%
% Builds a joint unsupervised UMAP embedding of all DeePhys units plus the
% externally labelled units, then uses graph label propagation to assign
% cell-type labels to the DeePhys units. External unit labels serve as seeds.
%
% This mirrors the standard classifyUnits pipeline, but swaps the Louvain-
% derived training labels for externally supplied ground truth (e.g. from
% patch-clamp or optogenetic identification).
%
% Normalization:
%   DeePhys features — per-group z-score (NormalizationVar, e.g. ChipID)
%                      via buildNormalizedFeatures, then joint global z-score.
%   External features — no per-group z-score (no metadata available);
%                       joint global z-score applied directly.
%   Both sides receive the same max(abs) scaling and feature group weighting.
%
% INPUTS:
%   ctc       - CellTypeClassifier with UnitDataArray set
%   wf_train  - (N_samples × N_train) raw unnormalised waveforms at sr_train
%   acg_train - (N_bins × N_train)   ACGs computed with the same ACGBinSize /
%               ACGLag as ctc.Parameters.Harmonization
%   y_train   - (1 × N_train) class labels (1=excitatory, 2=inhibitory)
%   sr_train  - waveform sampling rate of the external training data in Hz
%
% OPTIONAL NAME-VALUE:
%   ACGBinSize - bin size (s) used to compute acg_train. When provided,
%                it is checked against ctc.Parameters.Harmonization.ACGBinSize.
%   ACGLag            - half-width (s) used to compute acg_train. Checked against
%                       ctc.Parameters.Harmonization.ACGLag.
%   WaveformPreTrough - ms of waveform captured before the trough. Checked against
%                       ctc.Parameters.Harmonization.WaveformPreTrough.
%   WaveformPostTrough - ms of waveform captured after the trough. Checked against
%                        ctc.Parameters.Harmonization.WaveformPostTrough.
%
%   Omitting ACG / waveform provenance parameters issues a warning.
%   Use CellTypeClassifier.computeACGs() to generate compatible ACGs from
%   spike times.
%
% Sets:  ctc.UnitLabels            (1 = excitatory, 2 = inhibitory)
%        ctc.UnitConfidence        class probability after propagation
%        ctc.UnitGraphConnectivity number of direct external-unit graph neighbors
%        ctc.Reduction.ExternalTrain.DeePhys   (n_unique × NDims) DeePhys embedding
%        ctc.Reduction.ExternalTrain.External  (N_train  × NDims) external embedding
%        ctc.Reduction.ExternalTrain.ExternalLabels  y_train (for diagnostics)

arguments
    ctc       CellTypeClassifier
    wf_train  (:,:) double
    acg_train (:,:) double
    y_train   (1,:) double
    sr_train  (1,1) double
    opts.ACGBinSize         double = []
    opts.ACGLag             double = []
    opts.WaveformPreTrough  double = []
    opts.WaveformPostTrough double = []
end

assert(~isempty(ctc.UnitDataArray), ...
    'ctc.UnitDataArray is empty — set UnitDataArray before calling classifyUnitsWithExternalTrain().');

% ── ACG parameter validation ───────────────────────────────────────────────────
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
        ['Waveform window provenance unknown — verify wf_train was captured with at least ' ...
         '%.3f ms before and %.3f ms after the trough.'], ph.WaveformPreTrough, ph.WaveformPostTrough);
end

% ── Build DeePhys feature matrix (per-group z-scored, unique units only) ───────
ctc.buildNormalizedFeatures();
nf         = ctc.NormalizedFeatures;
X_dphy     = nf.X_pergroup;    % (n_unique × F) per-group z-scored
feat_names = nf.feat_names;    % (1 × F)
n_dphy     = size(X_dphy, 1);

% ── Build external train features ──────────────────────────────────────────────
[X_ext, feat_names_ext] = buildFeatureMatrix(ctc, wf_train, acg_train, sr_train);
n_ext = size(X_ext, 1);

assert(isequal(feat_names, feat_names_ext), ...
    ['Feature mismatch: external train (%d features) vs DeePhys (%d features).\n' ...
     'Ensure external ACGs use the same bin size / lag as ctc.Parameters.Harmonization.'], ...
    numel(feat_names_ext), numel(feat_names));

assert(n_ext == numel(y_train), ...
    'y_train length (%d) does not match number of external training units (%d).', ...
    numel(y_train), n_ext);

X_ext(isnan(X_ext)) = 0;

% ── Joint normalization ─────────────────────────────────────────────────────────
% DeePhys already has per-group z-score; external has none (no chip metadata).
% Global z-score is fit on the combined population so both sources land in the
% same feature space for the UMAP geometry to be meaningful.
X_combined = [X_dphy; X_ext];

[X_combined, ~, ~] = normalize(X_combined, 1, 'zscore');

% Remove columns that are NaN or Inf in any row
nan_cols = any(isnan(X_combined) | isinf(X_combined), 1);
X_combined(:, nan_cols) = [];
feat_names_trimmed = feat_names(~nan_cols);

% Scale by max(abs) over the full combined set
scale = max(abs(X_combined), [], 1);
scale(scale == 0) = 1;
X_combined = X_combined ./ scale;

% ── Dimension-normalised feature group weighting ───────────────────────────────
p_umap = ctc.Parameters.UMAP;
is_acg = startsWith(feat_names_trimmed, "FullACG") | startsWith(feat_names_trimmed, "ACG");
is_wf  = startsWith(feat_names_trimmed, "Waveform") | startsWith(feat_names_trimmed, "ReferenceWaveform");
n_acg  = sum(is_acg);
n_wf_f = sum(is_wf);
if n_acg > 0
    X_combined(:, is_acg) = X_combined(:, is_acg) * sqrt(p_umap.ACGWeight / n_acg);
end
if n_wf_f > 0
    X_combined(:, is_wf) = X_combined(:, is_wf) * sqrt(p_umap.WaveformWeight / n_wf_f);
end

% ── Unsupervised UMAP on combined DeePhys + external set ───────────────────────
N_all = n_dphy + n_ext;
if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_all)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_all);
else
    n_neighbors = p_umap.NNeighbors;
end

rng(ctc.Parameters.RNGSeed, 'twister');
fprintf('classifyUnitsWithExternalTrain: UMAP on %d units (%d DeePhys + %d external)\n', ...
    N_all, n_dphy, n_ext);

[reduction, umap_model, ~, ~] = run_umap(X_combined, ...
    'n_components', p_umap.NDims, ...
    'n_neighbors',  n_neighbors, ...
    'min_dist',     p_umap.MinDist, ...
    'spread',       p_umap.Spread, ...
    'method',       'Java', ...
    'verbose',      'none');

% ── Build row-stochastic graph from UMAP COO edge list ─────────────────────────
G = sparse(umap_model.head, umap_model.tail, umap_model.weights, N_all, N_all);
G = (G + G') / 2;
row_sums = sum(G, 2);
row_sums(row_sums == 0) = 1;
W = G ./ row_sums;

% ── Label propagation: external units are seeds, DeePhys units receive labels ──
p_conf = ctc.Parameters.Classification;

% Initialise: uninformed prior for DeePhys, clamped labels for external
F_label = ones(N_all, 2) * 0.5;
for i = 1:n_ext
    gi = n_dphy + i;
    if y_train(i) == 1
        F_label(gi, :) = [1, 0];
    else
        F_label(gi, :) = [0, 1];
    end
end

max_iter = p_conf.GraphMaxIter;
tol      = p_conf.GraphConvergenceTol;
for iter = 1:max_iter
    F_new = W * F_label;
    F_row_sum = sum(F_new, 2);
    F_row_sum(F_row_sum == 0) = 1;
    F_new = F_new ./ F_row_sum;
    % Clamp external seed labels
    for i = 1:n_ext
        gi = n_dphy + i;
        if y_train(i) == 1
            F_new(gi, :) = [1, 0];
        else
            F_new(gi, :) = [0, 1];
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

% ── Extract DeePhys labels, confidence, and graph connectivity ─────────────────
unique_labels     = nan(1, n_dphy);
unique_confidence = nan(1, n_dphy);
for i = 1:n_dphy
    p_exc = F_label(i, 1);
    p_inh = F_label(i, 2);
    if p_exc >= p_inh
        unique_labels(i)     = 1;
        unique_confidence(i) = p_exc;
    else
        unique_labels(i)     = 2;
        unique_confidence(i) = p_inh;
    end
end

% Number of direct external-unit neighbors in the UMAP graph for each DeePhys unit
graph_train_neigh = full(sum(G(1:n_dphy, n_dphy+1:end) > 0, 2))';  % (1 × n_dphy)

low_conn = sum(graph_train_neigh < 3);
if low_conn > 0.1 * n_dphy
    warning('CellTypeClassifier:sparseGraphNeighborhood', ...
        ['%d/%d DeePhys units have fewer than 3 external-unit neighbors in the ' ...
         'UMAP graph. Their confidence scores may reflect sparse connectivity ' ...
         'rather than robust classification.'], low_conn, n_dphy);
end

% Optional confidence threshold
if p_conf.UseConfidenceThreshold
    low_conf = unique_confidence < p_conf.ConfidenceThreshold;
    unique_labels(low_conf) = NaN;
    fprintf('Confidence threshold %.2f: %d units set to NaN\n', ...
        p_conf.ConfidenceThreshold, sum(low_conf));
end

% Broadcast unique-unit results back to all (unit × recording) rows
all_to_unique             = nf.all_to_unique;
ctc.UnitLabels            = unique_labels(all_to_unique);
ctc.UnitConfidence        = unique_confidence(all_to_unique);
ctc.UnitGraphConnectivity = graph_train_neigh(all_to_unique);

% Store combined embedding for diagnostics
ctc.Reduction.ExternalTrain = struct( ...
    'DeePhys',        reduction(1:n_dphy, :), ...
    'External',       reduction(n_dphy+1:end, :), ...
    'ExternalLabels', y_train);

fprintf('Classified %d excitatory, %d inhibitory DeePhys units (external train, n=%d)\n', ...
    sum(ctc.UnitLabels == 1, 'omitnan'), sum(ctc.UnitLabels == 2, 'omitnan'), n_ext);
end
