function buildNormalizedFeatures(ctc)
% BUILDNORMALIZEDFEATURES  Extract features and apply per-group z-score; cache result.
%
% Shared preprocessing step called by both generateTrainLabels and classifyUnits.
% On second call returns immediately — idempotent by design. Construct a new
% CellTypeClassifier instance if you need to re-run with different parameters.
%
% What this does (once):
%   1. Deduplicate UnitDataArray to one representative per unique UnitID
%   2. Extract harmonized waveforms + ACGs (three-level cache via getOrExtract)
%   3. Build feature matrix (waveform + ACG features, common grid)
%   4. Fill NaN → 0
%   5. Per-group z-score within each NormalizationVar group
%   6. Compute training-culture subset mask (TrainingCultureIdx)
%
% Does NOT perform global z-score, NaN-column removal, or scaling — those are
% done in generateTrainLabels (computeGlobalNormalization) and stored in
% NormalizationParams for classifyUnits / external units. They are fit on ALL
% unique units (transductive; label-free, so this does not leak labels), not on
% the training-culture subset.
%
% Sets ctc.NormalizedFeatures (see classdef for field descriptions).
% Sets ctc.HarmonizedWaveforms, ctc.HarmonizedACGs, ctc.HarmonizedSR.

arguments
    ctc CellTypeClassifier
end

if ~isempty(ctc.NormalizedFeatures)
    return  % already computed for this instance
end

p_umap = ctc.Parameters.UMAP;

% -- Deduplicate: one representative UnitData per unique UnitID ---------------
[unique_ud, all_to_unique, ~] = ctc.uniqueUnitMap();
n_unique     = numel(unique_ud);
n_units_full = numel(ctc.UnitDataArray);

unique_to_rep = zeros(1, n_unique);
for u = 1:n_unique
    rows = find(all_to_unique == u);
    unique_to_rep(u) = rows(1);
end

% -- Extract harmonized waveforms + ACGs from unique units only ---------------
[wf, acg, sr] = ctc.getOrExtract(unique_ud);

% -- ACG denoising (Harmonization.ACGSmoothing: "auto" = split-half choice, or width in ms) --
ph = ctc.Parameters.Harmonization;
sm = 0;
if isfield(ph, 'ACGSmoothing'), sm = ph.ACGSmoothing; end
sm_info = struct();
if (isstring(sm) || ischar(sm)) && lower(string(sm)) == "auto"
    [sigma_ms, sm_info] = selectACGSmoothing({unique_ud.SpikeTimes}, acg, ph.ACGBinSize, ph.ACGLag);
    fprintf('ACG smoothing: sigma %.2f ms (split-half optimum %.2f ms on %d units, scaled x%.2f to the feature ACGs)\n', ...
        sigma_ms, sm_info.sigma_half_ms, sm_info.n_units, sm_info.scale);
else
    sigma_ms = double(sm);
end
acg = smoothACG(acg, sigma_ms / 1000 / ph.ACGBinSize);

[X_raw, feat_names, aligned_wf, norm_acgs] = buildFeatureMatrix(ctc, wf, acg, sr);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf(:, all_to_unique);
ctc.HarmonizedACGs      = norm_acgs(:, all_to_unique);
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% -- Fill NaN before per-group z-score (missing features treated as mean) -----
X_raw(isnan(X_raw)) = 0;

% -- Map unique units to FeatureStore rows (order-independent) ----------------
% unique_to_rep contains UnitDataArray indices, but FeatureStore columns
% (NormalizationVar, culture mask) are in FeatureStore order. Build a mapping
% from each unique unit to its FeatureStore row via UnitID+RecordingID matching.
ud = ctc.UnitDataArray;
rep_ud_ids  = string({ud(unique_to_rep).UnitID});
rep_rec_ids = string({ud(unique_to_rep).RecordingID});
rep_keys    = rep_ud_ids + "|" + rep_rec_ids;
fs_keys     = string(ctc.FeatureStore.UnitTable.UnitID(:))' + "|" + ...
              string(ctc.FeatureStore.UnitTable.RecordingID(:))';
[~, unique_to_fs] = ismember(rep_keys, fs_keys);  % (1 x n_unique) row vector

% -- Per-group z-score on all unique units ------------------------------------
norm_var = string(p_umap.NormalizationVar);
if isempty(norm_var) || strlength(norm_var) == 0
    norm_var = "";   % no per-group z-score requested
end
% Optional bridge reference {field, value}: each group is z-scored with the mean / SD of its
% reference units only (e.g. wells of one composition present on every plate), so groups that
% differ in composition keep that difference while technical offsets are removed.
ref = {};
if isfield(p_umap, 'NormalizationReference'), ref = p_umap.NormalizationReference; end
group_stats = struct([]); group_col_unique = [];
if norm_var ~= "" && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col        = ctc.FeatureStore.UnitTable.(norm_var);
    group_col_unique = group_col(unique_to_fs);
    is_ref = true(numel(group_col_unique), 1);
    if ~isempty(ref)
        ref_col = ctc.FeatureStore.UnitTable.(ref{1});
        is_ref  = string(ref_col(unique_to_fs(:))) == string(ref{2});
    end
    % ConstantAs NaN as MATLAB's normalize: a feature without variance in any group is
    % dropped by computeGlobalNormalization
    [X_raw, fallback, group_stats] = zscoreByGroup(X_raw, group_col_unique, is_ref, 20, ConstantAs=NaN);
    if ~isempty(fallback)
        warning('CellTypeClassifier:fewReferenceUnits', ['%d %s group(s) have fewer than 20 reference ' ...
            'units (%s = %s) and are z-scored with their own statistics.'], numel(fallback), norm_var, ...
            ref{1}, string(ref{2}));
    end
elseif norm_var ~= ""
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in UnitTable — per-group z-score skipped.', norm_var);
end

% -- Training culture subset mask ---------------------------------------------
if ~isempty(p_umap.TrainingCultureIdx)
    subset_mask_full = CellTypeClassifier.buildCultureSubsetMask( ...
        ctc.FeatureStore, p_umap.TrainingCultureIdx, ctc.Parameters.CultureKeys);
    subset_mask = subset_mask_full(unique_to_fs);
else
    subset_mask = true(1, n_unique);
end

% -- Cache result -------------------------------------------------------------
nf               = struct();
nf.X_pergroup    = X_raw;
nf.feat_names    = feat_names;
nf.unique_ud     = unique_ud;
nf.all_to_unique = all_to_unique;
nf.unique_to_rep = unique_to_rep;
nf.unique_to_fs  = unique_to_fs;   % (1 x n_unique) FeatureStore row index for each unique unit
nf.subset_mask   = subset_mask;
nf.acg_sigma_ms  = sigma_ms;      % ACG smoothing applied (ms); external units get the same
nf.acg_smoothing = sm_info;
nf.n_unique      = n_unique;
nf.n_units_full  = n_units_full;
nf.group_stats   = group_stats;        % per-group mean / SD of the z-score (empty: none applied)
nf.group_values  = group_col_unique;   % group of each unique unit, to standardise new feature rows
ctc.NormalizedFeatures = nf;
end
