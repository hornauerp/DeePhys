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
% fit on the training subset in generateTrainLabels and stored in
% NormalizationParams. classifyUnits then applies those stored parameters to
% all units. This split is intentional: the global normalization must be fit
% on the training subset only.
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
[X_raw, feat_names, aligned_wf, norm_acgs] = buildFeatureMatrix(ctc, wf, acg, sr);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf;
ctc.HarmonizedACGs      = norm_acgs;
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% -- Fill NaN before per-group z-score (missing features treated as mean) -----
X_raw(isnan(X_raw)) = 0;

% -- Per-group z-score on all unique units ------------------------------------
norm_var = p_umap.NormalizationVar;
if ~isempty(norm_var) && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col        = ctc.FeatureStore.UnitTable.(norm_var);
    group_col_unique = group_col(unique_to_rep);
    [~, ~, iG]       = unique(string(group_col_unique), 'stable');
    for g = 1:max(iG)
        mask = (iG == g);
        X_raw(mask, :) = normalize(X_raw(mask, :));
    end
elseif ~isempty(norm_var)
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in UnitTable — per-group z-score skipped.', norm_var);
end

% -- Training culture subset mask ---------------------------------------------
if ~isempty(p_umap.TrainingCultureIdx)
    subset_mask_full = CellTypeClassifier.buildCultureSubsetMask( ...
        ctc.FeatureStore, p_umap.TrainingCultureIdx, ctc.Parameters.CultureKeys);
    subset_mask = subset_mask_full(unique_to_rep);
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
nf.subset_mask   = subset_mask;
nf.n_unique      = n_unique;
nf.n_units_full  = n_units_full;
ctc.NormalizedFeatures = nf;
end
