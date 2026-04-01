function identifyResponsiveUnits(ctc, metadata_filter)
% IDENTIFYRESPONSIVEUNITS  Bootstrap test to identify units with significant firing rate changes.
%
% Groups units by culture (identity keys in MetadataTable), then for each
% culture compares pre/post spike rates with a bootstrap permutation test.
% Units showing a significant change in the configured direction are flagged
% as responsive candidates.
%
% INPUTS:
%   ctc             - CellTypeClassifier (new API with FeatureStore + UnitDataArray)
%   metadata_filter - (optional) {field, value} cell pair to restrict which cultures
%                     contribute responsive units (e.g. {'AAV', 128})
%
% Direction is controlled by Parameters.Bootstrap.Direction:
%   "increase" - flag units with significant firing rate increase (default)
%                → inhibitory candidates in standard DREADD assay designs
%   "decrease" - flag units with significant firing rate decrease
%                → excitatory candidates (e.g. silencing excitatory population)
%   "both"     - flag units responding significantly in either direction
%
% Note: generateTrainLabels() treats all responsive units as inhibitory
% candidates. When using Direction="decrease", set TrainLabels manually
% (see run_clf_metadata.m) to reflect the correct biological assignment.
%
% Sets: ctc.ResponsiveUnitIdx, ctc.ResponsiveUnitDirection

arguments
    ctc             CellTypeClassifier
    metadata_filter cell = {}
end

% ── Validate metadata_filter ─────────────────────────────────────────────────
if ~isempty(metadata_filter)
    assert(iscell(metadata_filter) && numel(metadata_filter) == 2, ...
        'metadata_filter must be a 2-element cell: {fieldName, value}. Got %d elements.', ...
        numel(metadata_filter));
    filter_field = metadata_filter{1};
    filter_val   = metadata_filter{2};
else
    filter_field = '';
    filter_val   = [];
end

p   = ctc.Parameters.Bootstrap;
fs  = ctc.FeatureStore;
ud  = ctc.UnitDataArray;
direction = lower(string(p.Direction));

% Build culture ID per unit (from MetadataTable via RecordingID join)
meta            = fs.MetadataTable;
unit_rec_ids    = fs.UnitTable.RecordingID;
culture_ids     = FeatureStore.getCultureIDsForUnits( ...
    unit_rec_ids, meta, ctc.Parameters.CultureKeys);
unique_cultures = unique(culture_ids, 'stable');

% Match UnitDataArray to FeatureStore order by UnitID
unit_ids_in_table = fs.UnitTable.UnitID;
ud_ids            = string({ud.UnitID})';
[~, ud_order]     = ismember(unit_ids_in_table, ud_ids);

N_units               = numel(unit_ids_in_table);
responsive_idx        = false(1, N_units);
responsive_direction  = repmat("none", 1, N_units);

% Validate filter field exists in metadata (only once, outside culture loop)
if ~isempty(filter_field)
    assert(ismember(filter_field, meta.Properties.VariableNames), ...
        'metadata_filter field "%s" not found in MetadataTable. Available fields: %s', ...
        filter_field, strjoin(meta.Properties.VariableNames, ', '));
end

n_skipped = 0;

for c = 1:numel(unique_cultures)
    cid        = unique_cultures(c);
    unit_mask  = culture_ids == cid;
    unit_rows  = find(unit_mask);

    if isempty(unit_rows)
        continue
    end

    % Apply optional metadata filter
    if ~isempty(filter_field)
        rec_id_mask = ismember(meta.RecordingID, unit_rec_ids(unit_mask));
        if ~any(rec_id_mask)
            n_skipped = n_skipped + 1;
            continue
        end
        first_rec = meta(find(rec_id_mask, 1), :);
        rec_val   = first_rec.(filter_field);
        % Normalise to comparable types before equality test
        if isnumeric(filter_val) && isnumeric(rec_val)
            match = isequal(rec_val, filter_val);
        else
            match = isequal(string(rec_val), string(filter_val));
        end
        if ~match
            n_skipped = n_skipped + 1;
            continue
        end
    end

    % Build binned spike matrix (N_units × N_bins) for pre and post windows
    ud_indices = ud_order(unit_rows);
    valid      = ud_indices > 0;
    ud_indices = ud_indices(valid);
    unit_rows  = unit_rows(valid);
    if isempty(ud_indices), continue, end

    pre_mat  = CellTypeClassifier.binnedSpikeMatrix( ...
        ud(ud_indices), p.PreCutout, p.BinSize);
    post_mat = CellTypeClassifier.binnedSpikeMatrix( ...
        ud(ud_indices), p.PostCutout, p.BinSize);

    response = bootstrapFiringRateResponse(pre_mat, post_mat, p.NIter, p.Alpha);

    switch direction
        case "increase"
            flag_rows = unit_rows(response.increase);
            responsive_idx(flag_rows)       = true;
            responsive_direction(flag_rows) = "increase";
        case "decrease"
            flag_rows = unit_rows(response.decrease);
            responsive_idx(flag_rows)       = true;
            responsive_direction(flag_rows) = "decrease";
        case "both"
            inc_rows = unit_rows(response.increase);
            dec_rows = unit_rows(response.decrease);
            responsive_idx(inc_rows)        = true;
            responsive_idx(dec_rows)        = true;
            responsive_direction(inc_rows)  = "increase";
            responsive_direction(dec_rows)  = "decrease";
        otherwise
            error('CellTypeClassifier:invalidDirection', ...
                'Parameters.Bootstrap.Direction must be "increase", "decrease", or "both". Got "%s".', ...
                p.Direction);
    end
end

ctc.ResponsiveUnitIdx       = responsive_idx;
ctc.ResponsiveUnitDirection = responsive_direction;

if n_skipped > 0
    fprintf('Responsive units (%s): %i / %i (%.1f%%)  [%i / %i cultures skipped by filter]\n', ...
        p.Direction, sum(responsive_idx), N_units, 100 * mean(responsive_idx), ...
        n_skipped, numel(unique_cultures));
else
    fprintf('Responsive units (%s): %i / %i (%.1f%%)\n', ...
        p.Direction, sum(responsive_idx), N_units, 100 * mean(responsive_idx));
end
end
