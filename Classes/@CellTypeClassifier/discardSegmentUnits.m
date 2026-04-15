function discardSegmentUnits(ctc)
% DISCARDSEGMENTUNITS  Remove non-baseline recording segments from the workflow.
%
% In dose-response experiments, the same neuron appears in multiple recordings
% (segments at different doses). identifyResponsiveUnits() uses all segments
% to establish the dose-response curve, but subsequent steps (generateTrainLabels,
% classifyUnits) should operate on baseline (representative) units only.
%
% This method:
%   1. Stores the full UnitDataArray in OriginalUnitDataArray (for FullACG
%      recomputation from combined parent spike trains)
%   2. Identifies baseline rows (GroupingVar == GroupingValues)
%   3. Filters UnitDataArray, FeatureStore.UnitTable, and ResponsiveUnit* arrays
%   4. Clears cached state that depends on array dimensions
%
% Call AFTER identifyResponsiveUnits() and BEFORE generateTrainLabels().
%
% NOTE: This modifies FeatureStore.UnitTable in-place (FeatureStore is a handle).
%       MetadataTable and RecordingTable are preserved unfiltered so that
%       dose-response diagnostics can still access all recordings.

arguments
    ctc CellTypeClassifier
end

% Deep-copy FeatureStore to avoid mutating a shared handle.
% FeatureStore is a handle class; in-place mutation would be visible to any
% other CellTypeClassifier instance (or caller) holding the same reference.
% FeatureStore.subset() creates a new instance and copies MetadataFields properly.
all_rec_ids = string(ctc.FeatureStore.MetadataTable.RecordingID);
ctc.FeatureStore = ctc.FeatureStore.subset(all_rec_ids);

% Store original for FullACG recomputation via computeParentACGs
ctc.OriginalUnitDataArray = ctc.UnitDataArray;

group_field = string(ctc.Parameters.UMAP.GroupingVar);
baseline_gv = ctc.Parameters.UMAP.GroupingValues;
meta        = ctc.FeatureStore.MetadataTable;
ut          = ctc.FeatureStore.UnitTable;
n_rows      = height(ut);

if ~ismember(group_field, meta.Properties.VariableNames)
    fprintf('discardSegmentUnits: GroupingVar "%s" not in MetadataTable — no segments to discard.\n', group_field);
    return
end

% -- Determine which recordings are baseline ----------------------------------
rec_ids     = string(ut.RecordingID);
unique_recs = unique(rec_ids, 'stable');
n_recs      = numel(unique_recs);
rec_is_baseline = false(n_recs, 1);

for ri = 1:n_recs
    rec_row = meta(meta.RecordingID == unique_recs(ri), :);
    if isempty(rec_row)
        rec_is_baseline(ri) = true;   % keep units with no metadata
        continue
    end
    gv = rec_row.(group_field)(1);
    if isnumeric(gv) && isnumeric(baseline_gv)
        rec_is_baseline(ri) = any(gv == baseline_gv);
    else
        rec_is_baseline(ri) = any(string(gv) == string(baseline_gv));
    end
end

% Map recording-level baseline flag to unit-level
[~, rec_idx] = ismember(rec_ids, unique_recs);
is_baseline  = rec_is_baseline(rec_idx)';   % (1 x n_rows) logical

n_after = sum(is_baseline);
if n_after == n_rows
    fprintf('discardSegmentUnits: all %d units are baseline — nothing to discard.\n', n_rows);
    return
end

% -- Filter UnitDataArray and FeatureStore.UnitTable ---------------------------
ctc.UnitDataArray          = ctc.UnitDataArray(is_baseline);
ctc.FeatureStore.UnitTable = ut(is_baseline, :);

% -- Filter ResponsiveUnit* arrays (FeatureStore order, same as UnitTable) -----
if ~isempty(ctc.ResponsiveUnitIdx)
    ctc.ResponsiveUnitIdx = ctc.ResponsiveUnitIdx(is_baseline);
end
if ~isempty(ctc.ResponsiveUnitDirection)
    ctc.ResponsiveUnitDirection = ctc.ResponsiveUnitDirection(is_baseline);
end
if ~isempty(ctc.ResponsiveStrength)
    ctc.ResponsiveStrength = ctc.ResponsiveStrength(is_baseline);
end
if ~isempty(ctc.CounterexampleUnitIdx)
    ctc.CounterexampleUnitIdx = ctc.CounterexampleUnitIdx(is_baseline);
end

% -- Clear cached state (must be recomputed on filtered arrays) ----------------
ctc.NormalizedFeatures  = [];
ctc.CachedExtraction    = [];
ctc.HarmonizedWaveforms = [];
ctc.HarmonizedACGs      = [];

n_unique_before = numel(unique(string({ctc.OriginalUnitDataArray.UnitID})));
n_unique_after  = numel(unique(string({ctc.UnitDataArray.UnitID})));
fprintf('discardSegmentUnits: %d -> %d rows (%d -> %d unique units, %d segments discarded)\n', ...
    n_rows, n_after, n_unique_before, n_unique_after, n_rows - n_after);
end
