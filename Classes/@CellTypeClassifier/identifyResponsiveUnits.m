function identifyResponsiveUnits(ctc, metadata_filter)
% IDENTIFYRESPONSIVEUNITS  Identify units with significant firing rate changes.
%
% Two methods controlled by Parameters.Bootstrap.GroundTruthMethod:
%
%   "two_window" (default)
%       Bootstrap permutation test comparing binned spike counts in a pre
%       window vs a post window within a single recording.
%       Parameters used: PreCutout, PostCutout, BinSize, NIter, Alpha.
%
%   "full_curve"
%       Per-unit Spearman rank correlation between firing rate and the
%       GroupingVar (e.g. Concentration) across all recordings in the
%       culture. The same UnitID must be present across recordings (i.e.
%       cross-recording unit tracking is required).
%       Parameters used: FullCurveAlpha, MinRecordings.
%       PreCutout / PostCutout / BinSize / NIter / Alpha are NOT used.
%       Cultures with fewer than MinRecordings recordings fall back to
%       "two_window" with a printed notice.
%
% Direction is controlled by Parameters.Bootstrap.Direction:
%   "increase" - units whose FR increases with GroupingVar (default)
%   "decrease" - units whose FR decreases with GroupingVar
%   "both"     - either direction
%
% Optional FDR correction: Parameters.Bootstrap.UseFDR (Benjamini-Hochberg).
% Optional strength filter: Parameters.Bootstrap.MinResponsiveStrength.
%
% INPUTS:
%   ctc             - CellTypeClassifier
%   metadata_filter - (optional) {field, value} cell pair to restrict which
%                     cultures contribute responsive units

arguments
    ctc             CellTypeClassifier
    metadata_filter cell = {}
end

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

p           = ctc.Parameters.Bootstrap;
fs          = ctc.FeatureStore;
ud          = ctc.UnitDataArray;
direction   = lower(string(p.Direction));
method      = lower(string(p.GroundTruthMethod));
group_field = string(ctc.Parameters.UMAP.GroupingVar);

meta            = fs.MetadataTable;
unit_rec_ids    = fs.UnitTable.RecordingID;
culture_ids     = FeatureStore.getCultureIDsForUnits( ...
    unit_rec_ids, meta, ctc.Parameters.CultureKeys);
unique_cultures = unique(culture_ids, 'stable');

unit_ids_in_table = fs.UnitTable.UnitID;

% Map each FeatureStore row to the correct UnitData entry using
% (UnitID, RecordingID) as a compound key.
% UnitID alone is ambiguous when the same unit spans multiple recordings.
ud_keys    = string({ud.UnitID})' + "|" + string({ud.RecordingID})';
table_keys = string(unit_ids_in_table) + "|" + string(unit_rec_ids);
[~, ud_order] = ismember(table_keys, ud_keys);

N_units              = numel(unit_ids_in_table);
responsive_idx       = false(1, N_units);
responsive_direction = repmat("none", 1, N_units);
responsive_strength  = zeros(1, N_units);
all_pvals            = nan(1, N_units);

if ~isempty(filter_field)
    assert(ismember(filter_field, meta.Properties.VariableNames), ...
        'metadata_filter field "%s" not found in MetadataTable. Available: %s', ...
        filter_field, strjoin(meta.Properties.VariableNames, ', '));
end

n_skipped = 0;

for c = 1:numel(unique_cultures)
    cid       = unique_cultures(c);
    unit_mask = culture_ids == cid;
    unit_rows = find(unit_mask);
    if isempty(unit_rows), continue, end

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

    ud_indices = ud_order(unit_rows);
    valid      = ud_indices > 0;
    ud_indices = ud_indices(valid);
    unit_rows  = unit_rows(valid);
    if isempty(ud_indices), continue, end

    rec_ids_here   = unique(unit_rec_ids(unit_mask), 'stable');
    use_full_curve = strcmp(method, 'full_curve') && ...
                     numel(rec_ids_here) >= p.MinRecordings && ...
                     ismember(group_field, meta.Properties.VariableNames);

    if use_full_curve
        % -- Full curve: per-unit Spearman dose-response ----------------------
        [flag_inc, flag_dec, strengths, pvals] = fullCurveResponsive( ...
            ud, ud_indices, unit_rows, unit_ids_in_table, ...
            meta, rec_ids_here, unit_rec_ids, unit_mask, ud_order, ...
            group_field, p.FullCurveAlpha);

        all_pvals(unit_rows)           = pvals;
        responsive_strength(unit_rows) = strengths;

        switch direction
            case "increase"
                flag_rows = unit_rows(flag_inc);
                responsive_idx(flag_rows)       = true;
                responsive_direction(flag_rows) = "increase";
            case "decrease"
                flag_rows = unit_rows(flag_dec);
                responsive_idx(flag_rows)       = true;
                responsive_direction(flag_rows) = "decrease";
            case "both"
                responsive_idx(unit_rows(flag_inc))       = true;
                responsive_idx(unit_rows(flag_dec))       = true;
                responsive_direction(unit_rows(flag_inc)) = "increase";
                responsive_direction(unit_rows(flag_dec)) = "decrease";
            otherwise
                error('CellTypeClassifier:invalidDirection', ...
                    'Bootstrap.Direction must be "increase", "decrease", or "both". Got "%s".', ...
                    p.Direction);
        end

    else
        % -- Two window: bootstrap pre/post comparison ------------------------
        if strcmp(method, 'full_curve')
            fprintf('Culture %d: %d / %d recordings (< MinRecordings=%d) or GroupingVar "%s" missing — falling back to two_window\n', ...
                c, numel(rec_ids_here), p.MinRecordings, p.MinRecordings, group_field);
        end

        pre_mat  = CellTypeClassifier.binnedSpikeMatrix(ud(ud_indices), p.PreCutout,  p.BinSize);
        post_mat = CellTypeClassifier.binnedSpikeMatrix(ud(ud_indices), p.PostCutout, p.BinSize);
        response = bootstrapFiringRateResponse(pre_mat, post_mat, p.NIter, p.Alpha);

        if isfield(response, 'p_values') && ~isempty(response.p_values)
            all_pvals(unit_rows) = response.p_values;
        end
        responsive_strength(unit_rows) = response.strength;

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
end

ctc.ResponsiveUnitIdx       = responsive_idx;
ctc.ResponsiveUnitDirection = responsive_direction;
ctc.ResponsiveStrength      = responsive_strength;

% -- FDR correction across all unit p-values ---------------------------------
if p.UseFDR
    valid_pval = ~isnan(all_pvals);
    if sum(valid_pval) >= 2
        pvals_test             = all_pvals(valid_pval);
        test_idx               = find(valid_pval);
        [sorted_p, sort_order] = sort(pvals_test);
        m         = numel(sorted_p);
        bh_thresh = (1:m)' / m * p.FDRLevel;
        max_sig   = find(sorted_p(:) <= bh_thresh, 1, 'last');

        if ~isempty(max_sig)
            sig_order   = sort_order(1:max_sig);
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
            'UseFDR=true but no p_values available. Check method output.');
    end
end

% -- Strength-based filtering ------------------------------------------------
if p.MinResponsiveStrength > 0
    weak      = ctc.ResponsiveUnitIdx & (responsive_strength < p.MinResponsiveStrength);
    n_demoted = sum(weak);
    ctc.ResponsiveUnitIdx(weak) = false;
    if n_demoted > 0
        fprintf('Demoted %d weak responders (strength < %.2f)\n', n_demoted, p.MinResponsiveStrength);
    end
end

% Count unique units (same UnitID spans multiple rows in full_curve mode)
n_responsive = numel(unique(unit_ids_in_table(ctc.ResponsiveUnitIdx)));
n_total      = numel(unique(unit_ids_in_table));

if n_skipped > 0
    fprintf('Responsive units (%s): %i / %i (%.1f%%)  [%i / %i cultures skipped by filter]\n', ...
        p.Direction, n_responsive, n_total, 100 * n_responsive / n_total, ...
        n_skipped, numel(unique_cultures));
else
    fprintf('Responsive units (%s): %i / %i (%.1f%%)\n', ...
        p.Direction, n_responsive, n_total, 100 * n_responsive / n_total);
end
end

%% ── Helper: per-unit full-curve dose-response ─────────────────────────────

function [flag_inc, flag_dec, strengths, pvals] = fullCurveResponsive( ...
        ud, ud_indices, unit_rows, unit_ids_in_table, ...
        meta, rec_ids_here, unit_rec_ids, unit_mask, ud_order, ...
        group_field, alpha)
% For each unique UnitID in this culture, build a firing-rate vector across
% recordings (one FR per recording the unit appears in), then test its
% Spearman correlation against the GroupingVar.
%
% Vectorized: pre-extracts all RecordingIDs and FRs in one pass, builds a
% (n_unique_uids × n_valid_recs) FR matrix, then calls corr on the whole
% matrix at once. Only units with missing recordings fall back to a scalar loop.

    n_rows    = numel(unit_rows);
    flag_inc  = false(1, n_rows);
    flag_dec  = false(1, n_rows);
    strengths = zeros(1, n_rows);
    pvals     = ones(1, n_rows);

    % ── GroupingVar value per recording ──────────────────────────────────────
    n_recs     = numel(rec_ids_here);
    group_vals = nan(1, n_recs);
    for ri = 1:n_recs
        rec_row = meta(meta.RecordingID == rec_ids_here(ri), :);
        if isempty(rec_row); continue; end
        gv = rec_row.(group_field)(1);
        if isnumeric(gv)
            group_vals(ri) = gv;
        else
            gv_num = str2double(string(gv));
            if ~isnan(gv_num); group_vals(ri) = gv_num; end
        end
    end

    valid_rec_mask = ~isnan(group_vals);
    if sum(valid_rec_mask) < 2 || numel(unique(group_vals(valid_rec_mask))) < 2
        return
    end

    valid_rec_ids   = rec_ids_here(valid_rec_mask);
    gv_valid        = group_vals(valid_rec_mask)';   % (n_valid_recs × 1)
    n_valid_recs    = numel(valid_rec_ids);

    % ── Pre-extract RecordingID and FR for every UnitData in this culture ────
    % Single vectorized pass — no per-unit string conversion in a loop.
    all_rec_ids = string({ud(ud_indices).RecordingID})';          % (n_rows × 1)
    all_durs    = [ud(ud_indices).RecordingDuration]';             % (n_rows × 1)
    all_nspk    = cellfun(@numel, {ud(ud_indices).SpikeTimes})';  % (n_rows × 1)
    all_frs     = all_nspk ./ max(all_durs, eps);                 % (n_rows × 1)

    % ── Build FR matrix: (n_unique_uids × n_valid_recs) ──────────────────────
    row_unit_ids    = unit_ids_in_table(unit_rows);
    unique_uids     = unique(row_unit_ids, 'stable');
    n_unique        = numel(unique_uids);
    fr_matrix       = nan(n_unique, n_valid_recs);

    valid_rec_str = string(valid_rec_ids);
    for u = 1:n_unique
        local_mask  = row_unit_ids == unique_uids(u);  % rows of this unit in culture
        uid_frs     = all_frs(local_mask);
        uid_rec_ids = all_rec_ids(local_mask);
        for ri = 1:n_valid_recs
            hit = uid_rec_ids == valid_rec_str(ri);
            if any(hit)
                fr_matrix(u, ri) = uid_frs(find(hit, 1));
            end
        end
    end

    % ── Vectorized corr for fully-observed units ──────────────────────────────
    rhos   = nan(n_unique, 1);
    p_incs = ones(n_unique, 1);
    p_decs = ones(n_unique, 1);

    full_mask = all(~isnan(fr_matrix), 2);
    if any(full_mask)
        fr_full = fr_matrix(full_mask, :)';   % (n_valid_recs × n_full_units)
        [r,  pr] = corr(gv_valid, fr_full, 'Type', 'Spearman', 'Tail', 'right');
        [~,  pl] = corr(gv_valid, fr_full, 'Type', 'Spearman', 'Tail', 'left');
        rhos(full_mask)   = r';
        p_incs(full_mask) = pr';
        p_decs(full_mask) = pl';
    end

    % Scalar fallback for units missing some recordings
    partial_idx = find(~full_mask & sum(~isnan(fr_matrix), 2) >= 2);
    for u = partial_idx'
        valid_pts = ~isnan(fr_matrix(u, :));
        gv_fit    = gv_valid(valid_pts);
        fr_fit    = fr_matrix(u, valid_pts)';
        try
            [rhos(u),  p_incs(u)] = corr(gv_fit, fr_fit, 'Type', 'Spearman', 'Tail', 'right');
            [~,        p_decs(u)] = corr(gv_fit, fr_fit, 'Type', 'Spearman', 'Tail', 'left');
        catch
        end
    end

    % ── Map results back to unit_rows ─────────────────────────────────────────
    for u = 1:n_unique
        local_idx            = find(row_unit_ids == unique_uids(u));
        strengths(local_idx) = abs(rhos(u));
        pvals(local_idx)     = min(p_incs(u), p_decs(u));
        flag_inc(local_idx)  = p_incs(u) < alpha;
        flag_dec(local_idx)  = p_decs(u) < alpha;
    end
end
