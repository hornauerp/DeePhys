function identifyResponsiveUnits(ctc, metadata_filter)
% IDENTIFYRESPONSIVEUNITS  Bootstrap test to identify units with significant firing rate changes.
%
% Groups units by culture (identity keys in MetadataTable), then for each
% culture runs responsive unit detection. Two methods are available:
%
%   "two_window" (default) — Compares pre/post spike rates within a single recording
%       using a bootstrap permutation test. Units with a significant change in
%       the configured direction are flagged as responsive candidates.
%       Strength = |empirical difference| / null_std (effect size).
%
%   "full_curve" — For each culture, computes mean population firing rate per
%       recording and tests for monotonic dose-response using Spearman rank
%       correlation. Cultures with fewer than MinRecordings recordings fall back
%       to "two_window". Within responsive cultures, unit-level detection is
%       done via two_window. Strength = |culture_rho| * unit_effect_size.
%
% Direction is controlled by Parameters.Bootstrap.Direction:
%   "increase" - flag units with significant firing rate increase (default)
%   "decrease" - flag units with significant firing rate decrease
%   "both"     - flag units responding in either direction
%
% When UseFDR is true, a Benjamini-Hochberg correction is applied across all
% unit p-values (requires bootstrapFiringRateResponse to return .p_values).
% When MinResponsiveStrength > 0, units below this strength are demoted.
%
% INPUTS:
%   ctc             - CellTypeClassifier with FeatureStore + UnitDataArray
%   metadata_filter - (optional) {field, value} cell pair to restrict which cultures
%                     contribute responsive units (e.g. {'Concentration', 0})
%
% Sets: ctc.ResponsiveUnitIdx, ctc.ResponsiveUnitDirection, ctc.ResponsiveStrength

arguments
    ctc             CellTypeClassifier
    metadata_filter cell = {}
end

% -- Validate metadata_filter ------------------------------------------------
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

p         = ctc.Parameters.Bootstrap;
fs        = ctc.FeatureStore;
ud        = ctc.UnitDataArray;
direction = lower(string(p.Direction));
method    = lower(string(p.GroundTruthMethod));

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

N_units              = numel(unit_ids_in_table);
responsive_idx       = false(1, N_units);
responsive_direction = repmat("none", 1, N_units);
responsive_strength  = zeros(1, N_units);
all_pvals            = nan(1, N_units);   % for FDR correction

% Validate filter field (once, outside culture loop)
if ~isempty(filter_field)
    assert(ismember(filter_field, meta.Properties.VariableNames), ...
        'metadata_filter field "%s" not found in MetadataTable. Available: %s', ...
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

    % Match UnitDataArray indices to this culture's units
    ud_indices = ud_order(unit_rows);
    valid      = ud_indices > 0;
    ud_indices = ud_indices(valid);
    unit_rows  = unit_rows(valid);
    if isempty(ud_indices), continue, end

    % -- Full curve: culture-level dose-response prefilter -------------------
    rec_ids_here = unique(unit_rec_ids(unit_mask));
    use_full_curve = strcmp(method, 'full_curve') && numel(rec_ids_here) >= p.MinRecordings;

    culture_rho = NaN;
    if use_full_curve
        culture_rho = computeCultureDoseResponse(ud, ud_order, meta, ...
            rec_ids_here, unit_rec_ids, unit_mask, p);
        if isnan(culture_rho)
            fprintf('Culture %d: full_curve failed — falling back to two_window\n', c);
            use_full_curve = false;
        end
    end

    % -- Bootstrap test -------------------------------------------------------
    pre_mat  = CellTypeClassifier.binnedSpikeMatrix( ...
        ud(ud_indices), p.PreCutout, p.BinSize);
    post_mat = CellTypeClassifier.binnedSpikeMatrix( ...
        ud(ud_indices), p.PostCutout, p.BinSize);
    response = bootstrapFiringRateResponse(pre_mat, post_mat, p.NIter, p.Alpha);

    % Collect p-values for optional FDR correction
    if isfield(response, 'p_values') && ~isempty(response.p_values)
        for ui = 1:numel(unit_rows)
            all_pvals(unit_rows(ui)) = response.p_values(ui);
        end
    end

    % Assign unit-level strength (scale by culture rho for full_curve)
    for ui = 1:numel(unit_rows)
        if use_full_curve
            responsive_strength(unit_rows(ui)) = abs(culture_rho) * response.strength(ui);
        else
            responsive_strength(unit_rows(ui)) = response.strength(ui);
        end
    end

    % -- Direction-based responsive flagging ----------------------------------
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
                'Bootstrap.Direction must be "increase", "decrease", or "both". Got "%s".', ...
                p.Direction);
    end
end

ctc.ResponsiveUnitIdx       = responsive_idx;
ctc.ResponsiveUnitDirection = responsive_direction;
ctc.ResponsiveStrength      = responsive_strength;

% -- FDR correction across all unit p-values (Phase 2) ----------------------
if p.UseFDR
    valid_pval = ~isnan(all_pvals);
    if sum(valid_pval) >= 2
        pvals_test  = all_pvals(valid_pval);
        test_idx    = find(valid_pval);
        [sorted_p, sort_order] = sort(pvals_test);
        m = numel(sorted_p);
        bh_thresh = (1:m)' / m * p.FDRLevel;
        max_sig   = find(sorted_p(:) <= bh_thresh, 1, 'last');

        if ~isempty(max_sig)
            sig_order = sort_order(1:max_sig);
            significant = false(1, N_units);
            significant(test_idx(sig_order)) = true;
        else
            significant = false(1, N_units);
        end

        n_before = sum(responsive_idx);
        ctc.ResponsiveUnitIdx = responsive_idx & significant;
        fprintf('FDR correction (q=%.3f): %d -> %d responsive units\n', ...
            p.FDRLevel, n_before, sum(ctc.ResponsiveUnitIdx));
    else
        warning('CellTypeClassifier:fdrNoPValues', ...
            'UseFDR=true but no p_values returned from bootstrapFiringRateResponse. ' ...
            'Upgrade bootstrapFiringRateResponse or disable UseFDR.');
    end
end

% -- Optional strength-based filtering (Phase 1) ----------------------------
if p.MinResponsiveStrength > 0
    weak = ctc.ResponsiveUnitIdx & (responsive_strength < p.MinResponsiveStrength);
    n_demoted = sum(weak);
    ctc.ResponsiveUnitIdx(weak) = false;
    if n_demoted > 0
        fprintf('Demoted %d weak responders (strength < %.2f)\n', ...
            n_demoted, p.MinResponsiveStrength);
    end
end

if n_skipped > 0
    fprintf('Responsive units (%s): %i / %i (%.1f%%)  [%i / %i cultures skipped by filter]\n', ...
        p.Direction, sum(ctc.ResponsiveUnitIdx), N_units, ...
        100 * mean(ctc.ResponsiveUnitIdx), n_skipped, numel(unique_cultures));
else
    fprintf('Responsive units (%s): %i / %i (%.1f%%)\n', ...
        p.Direction, sum(ctc.ResponsiveUnitIdx), N_units, ...
        100 * mean(ctc.ResponsiveUnitIdx));
end
end

%% -- Helper: culture-level dose-response ------------------------------------

function rho = computeCultureDoseResponse(ud, ud_order, meta, rec_ids_here, ...
        unit_rec_ids, unit_mask, p)
% Compute Spearman rank correlation between population firing rate and
% GroupingVar (e.g. Concentration) across recordings within a culture.
% Returns NaN if computation fails or fewer than 2 unique GroupingVar values.

    group_field = 'Concentration';  % default; uses Parameters.UMAP.GroupingVar if available
    if ~ismember(group_field, meta.Properties.VariableNames)
        rho = NaN;
        return
    end

    n_recs     = numel(rec_ids_here);
    mean_rates = nan(1, n_recs);
    group_vals = nan(1, n_recs);

    for ri = 1:n_recs
        rid     = rec_ids_here(ri);
        rec_row = meta(meta.RecordingID == rid, :);
        if isempty(rec_row); continue; end
        group_vals(ri) = rec_row.(group_field)(1);

        % Units from this specific recording
        rec_unit_mask = unit_mask & (unit_rec_ids == rid);
        rec_rows      = find(rec_unit_mask);
        rec_ud_idx    = ud_order(rec_rows);
        rec_ud_idx    = rec_ud_idx(rec_ud_idx > 0);
        if isempty(rec_ud_idx); continue; end

        all_rates = zeros(1, numel(rec_ud_idx));
        for ui = 1:numel(rec_ud_idx)
            st = ud(rec_ud_idx(ui)).SpikeTimes;
            if isempty(st)
                all_rates(ui) = 0;
            else
                dur = max(max(st), p.PostCutout(2));
                all_rates(ui) = numel(st) / dur;
            end
        end
        mean_rates(ri) = mean(all_rates);
    end

    valid = ~isnan(mean_rates) & ~isnan(group_vals);
    if sum(valid) < 2 || numel(unique(group_vals(valid))) < 2
        rho = NaN;
        return
    end

    try
        rho = corr(group_vals(valid)', mean_rates(valid)', 'Type', 'Spearman');
    catch
        rho = NaN;
    end
end
