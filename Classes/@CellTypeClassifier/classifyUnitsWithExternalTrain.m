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
if isfield(nf, 'acg_sigma_ms') && nf.acg_sigma_ms > 0   % same ACG smoothing as the DeePhys units
    acg_train = smoothACG(acg_train, nf.acg_sigma_ms / 1000 / ph.ACGBinSize);
end
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
X_combined = applyFeatureGroupWeights(X_combined, feat_names_trimmed, p_umap);

% ── Unsupervised UMAP on combined DeePhys + external set ───────────────────────
N_all = n_dphy + n_ext;
if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_all)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_all);
else
    n_neighbors = p_umap.NNeighbors;
end

fprintf('classifyUnitsWithExternalTrain: joint graph on %d units (%d DeePhys + %d external)\n', ...
    N_all, n_dphy, n_ext);

% ── Joint unit graph (same construction as the main pipeline) ──────────────────
[G, reduction] = jointGraph(ctc, X_combined, n_neighbors, feat_names_trimmed);
if isempty(reduction), reduction = nan(N_all, 0); end

% ── Label propagation: external units are seeds, DeePhys units receive labels ──
% Same propagation as classifyUnits: units with no graph path to an external unit are
% left unlabelled instead of being called excitatory by the tie rule, and F is seed-count
% normalised (external label sets are usually unbalanced).
p_conf = ctc.Parameters.Classification;
[F_label, unreached_all, converged, iter] = CellTypeClassifier.propagateLabels(G, ...
    n_dphy + (1:n_ext), y_train, p_conf.GraphMaxIter, p_conf.GraphConvergenceTol);
if converged
    fprintf('Graph label propagation converged at iteration %d\n', iter);
else
    fprintf('Graph label propagation reached max iterations (%d)\n', p_conf.GraphMaxIter);
end
unreached = unreached_all(1:n_dphy);
if any(unreached)
    warning('CellTypeClassifier:unreached', ['%d units have no graph path to an ' ...
        'external training unit and are left unlabelled (NaN).'], sum(unreached));
end

% ── Extract DeePhys labels, confidence, and graph connectivity ─────────────────
unique_labels     = nan(1, n_dphy);
unique_confidence = nan(1, n_dphy);
for i = 1:n_dphy
    if unreached(i), continue, end      % label and confidence stay NaN
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
