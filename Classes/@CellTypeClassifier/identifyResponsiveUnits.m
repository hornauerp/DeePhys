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
%       Three-stage pipeline using the full dose-response curve:
%         Stage 1 — Monotonicity: |Spearman rho(dose rank, mean FR)| >= threshold
%         Stage 2 — Effect size:  fold change (baseline to max dose) >= threshold
%         Stage 3 — Bootstrap:    permutation test on baseline vs max dose
%       The same UnitID must be present across recordings (i.e.
%       cross-recording unit tracking is required).
%       Parameters used: MonotonicityThreshold, MinFoldChange,
%       ConfirmationAlpha, ConfirmationNIter, MinRecordings, BinSize.
%       Cultures with fewer than MinRecordings recordings fall back to
%       "two_window" with a printed notice.
%
%   "metadata"
%       Per-unit ground truth labels read directly from a column in
%       FeatureStore.UnitTable. No firing-rate tests are run. Both classes
%       have explicit ground truth (pure culture scenario).
%       Parameters used: LabelField, ResponsiveClassValue.
%         LabelField: column name in UnitTable (e.g. "CellType")
%         ResponsiveClassValue: value mapping to ResponsiveClassLabel
%           (e.g. "inhibitory"). Units with other non-empty values become
%           explicit counterexample ground truth (ctc.CounterexampleUnitIdx).
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

% -- Metadata method: read labels directly from UnitTable ----------------------
if strcmp(method, 'metadata')
    label_field = string(p.LabelField);
    resp_value  = string(p.ResponsiveClassValue);
    resp_label  = ctc.Parameters.TrainLabels.ResponsiveClassLabel;

    assert(~isempty(label_field) && label_field ~= "", ...
        'CellTypeClassifier:noLabelField', ...
        'Bootstrap.LabelField must be set when GroundTruthMethod is "metadata".');
    assert(~isempty(resp_value) && resp_value ~= "", ...
        'CellTypeClassifier:noResponsiveClassValue', ...
        'Bootstrap.ResponsiveClassValue must be set when GroundTruthMethod is "metadata".');
    assert(ismember(label_field, string(fs.UnitTable.Properties.VariableNames)), ...
        'CellTypeClassifier:labelFieldMissing', ...
        'LabelField "%s" not found in FeatureStore.UnitTable. Available: %s', ...
        label_field, strjoin(string(fs.UnitTable.Properties.VariableNames), ', '));

    raw_labels = string(fs.UnitTable.(label_field));
    N_total    = numel(raw_labels);

    resp_idx   = raw_labels == resp_value;
    % Counterexamples: any unit with a non-empty label that is NOT the responsive class
    ce_idx     = ~resp_idx & raw_labels ~= "" & raw_labels ~= "NaN";

    ctc.ResponsiveUnitIdx       = resp_idx(:)';
    ctc.ResponsiveUnitDirection = repmat("none", 1, N_total);
    ctc.ResponsiveStrength      = double(resp_idx(:)');  % 1 for labeled, 0 otherwise
    ctc.CounterexampleUnitIdx   = ce_idx(:)';

    n_resp = sum(resp_idx);
    n_ce   = sum(ce_idx);
    n_unlabeled = N_total - n_resp - n_ce;
    resp_class_name = 'inhibitory';
    if resp_label == 1; resp_class_name = 'excitatory'; end
    fprintf('Metadata labels: %d %s (class %d), %d counterexamples, %d unlabeled\n', ...
        n_resp, resp_class_name, resp_label, n_ce, n_unlabeled);
    return
end

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
            group_field, p);

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
        % -- Two window: compare FR between two recordings in the culture -----
        if strcmp(method, 'full_curve')
            fprintf('Culture %d: %d recordings < MinRecordings=%d or GroupingVar "%s" missing — falling back to two_window\n', ...
                c, numel(rec_ids_here), p.MinRecordings, group_field);
        end

        % Get GroupingVar value per recording in this culture
        n_recs_here = numel(rec_ids_here);
        gv_here     = nan(1, n_recs_here);
        if ismember(group_field, meta.Properties.VariableNames)
            for ri = 1:n_recs_here
                rec_row = meta(meta.RecordingID == rec_ids_here(ri), :);
                if isempty(rec_row); continue; end
                gv = rec_row.(group_field)(1);
                if isnumeric(gv)
                    gv_here(ri) = gv;
                else
                    gv_num = str2double(string(gv));
                    if ~isnan(gv_num); gv_here(ri) = gv_num; end
                end
            end
        end

        % Select pre and post recordings by GroupingVar value (or auto)
        valid_gv = ~isnan(gv_here);
        if ~isempty(p.PreGroupValue) && isnumeric(p.PreGroupValue)
            pre_ri = find(valid_gv & gv_here == p.PreGroupValue, 1);
        else
            [~, pre_ri] = min(gv_here);   % lowest GroupingVar = baseline
        end
        if ~isempty(p.PostGroupValue) && isnumeric(p.PostGroupValue)
            post_ri = find(valid_gv & gv_here == p.PostGroupValue, 1);
        else
            [~, post_ri] = max(gv_here);  % highest GroupingVar = max treatment
        end

        if isempty(pre_ri) || isempty(post_ri) || pre_ri == post_ri
            continue
        end

        pre_rec_id  = rec_ids_here(pre_ri);
        post_rec_id = rec_ids_here(post_ri);

        % Get matched pre/post UnitData for units present in both recordings
        tw_pre_idx  = [];  tw_post_idx  = [];  tw_rows_tw = [];
        row_unit_ids_c = unit_ids_in_table(unit_rows);
        unique_uids_c  = unique(row_unit_ids_c, 'stable');
        for u = 1:numel(unique_uids_c)
            uid       = unique_uids_c(u);
            uid_mask  = row_unit_ids_c == uid;
            uid_ud_raw    = ud_order(unit_rows(uid_mask));
            uid_rows_here = unit_rows(uid_mask);
            valid_ud      = uid_ud_raw > 0;
            uid_ud        = uid_ud_raw(valid_ud);
            uid_rows_here = uid_rows_here(valid_ud);  % filter in sync with uid_ud
            uid_rids  = string({ud(uid_ud).RecordingID})';
            pre_hit   = find(uid_rids == string(pre_rec_id),  1);
            post_hit  = find(uid_rids == string(post_rec_id), 1);
            if isempty(pre_hit) || isempty(post_hit); continue; end
            tw_pre_idx  = [tw_pre_idx,  uid_ud(pre_hit)];   %#ok<AGROW>
            tw_post_idx = [tw_post_idx, uid_ud(post_hit)];  %#ok<AGROW>
            % Use the post-recording row as the representative FeatureStore row
            post_row = uid_rows_here(post_hit);
            tw_rows_tw = [tw_rows_tw, post_row];             %#ok<AGROW>
        end

        if isempty(tw_pre_idx); continue; end

        % Build binned spike matrices using optional sub-window or full recording
        pre_dur  = ud(tw_pre_idx(1)).RecordingDuration;
        post_dur = ud(tw_post_idx(1)).RecordingDuration;
        pre_win  = p.PreCutout;
        post_win = p.PostCutout;
        if isempty(pre_win);  pre_win  = [0, pre_dur];  end
        if isempty(post_win); post_win = [0, post_dur]; end

        % Equalise number of bins so bootstrap permutation is valid
        n_bins   = floor(min(diff(pre_win), diff(post_win)) / p.BinSize);
        pre_win  = [pre_win(1),  pre_win(1)  + n_bins * p.BinSize];
        post_win = [post_win(1), post_win(1) + n_bins * p.BinSize];

        pre_mat  = CellTypeClassifier.binnedSpikeMatrix(ud(tw_pre_idx),  pre_win,  p.BinSize);
        post_mat = CellTypeClassifier.binnedSpikeMatrix(ud(tw_post_idx), post_win, p.BinSize);
        response = bootstrapFiringRateResponse(pre_mat, post_mat, p.NIter, p.Alpha);

        % Broadcast results to ALL rows of each tested unit (consistent
        % with full_curve, which flags every row of a responsive unit).
        % tw_rows_tw contains one representative row per unit; expand to
        % all rows sharing the same UnitID in this culture.
        tw_uids = unit_ids_in_table(tw_rows_tw);
        for ti = 1:numel(tw_rows_tw)
            all_uid_rows = unit_rows(row_unit_ids_c == tw_uids(ti));
            if isfield(response, 'p_values') && ~isempty(response.p_values)
                all_pvals(all_uid_rows) = response.p_values(ti);
            end
            responsive_strength(all_uid_rows) = response.strength(ti);
        end

        switch direction
            case "increase"
                flag_uids = tw_uids(response.increase);
            case "decrease"
                flag_uids = tw_uids(response.decrease);
            case "both"
                flag_uids = tw_uids(response.increase | response.decrease);
            otherwise
                error('CellTypeClassifier:invalidDirection', ...
                    'Bootstrap.Direction must be "increase", "decrease", or "both". Got "%s".', ...
                    p.Direction);
        end
        for fi = 1:numel(flag_uids)
            all_uid_rows = unit_rows(row_unit_ids_c == flag_uids(fi));
            responsive_idx(all_uid_rows) = true;
            if direction == "both"
                % Determine direction for this specific unit
                uid_ti = find(tw_uids == flag_uids(fi), 1);
                if response.increase(uid_ti)
                    responsive_direction(all_uid_rows) = "increase";
                else
                    responsive_direction(all_uid_rows) = "decrease";
                end
            else
                responsive_direction(all_uid_rows) = direction;
            end
        end
    end
end

ctc.ResponsiveUnitIdx       = responsive_idx;
ctc.ResponsiveUnitDirection = responsive_direction;
ctc.ResponsiveStrength      = responsive_strength;

% -- Store ResponsivenessDetail for diagnostics --------------------------------
% Build per-unit FR-vs-dose matrix from UnitData. For full_curve method this
% reuses the already-computed FR values via a lightweight pass over UnitData.
% For two_window, only pre and post FR are stored (2-point dose response).
ctc.ResponsivenessDetail = buildResponsivenessDetail(ctc, unit_ids_in_table, ...
    ud, ud_order, meta, group_field, method, p);

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

% -- Optional diagnostic ------------------------------------------------------
if ctc.Parameters.Diagnostics.Enable
    ctc.diagnosticResponsiveUnits();
end
end

%% ── Helper: build ResponsivenessDetail struct ─────────────────────────────

function detail = buildResponsivenessDetail(ctc, unit_ids_in_table, ud, ud_order, ...
        meta, group_field, method, p)
% Builds per-unit FR-vs-dose matrix for diagnostic plotting.
% For full_curve: extracts FR across all recordings sorted by GroupingVar.
% For two_window / metadata: stores 2-point FR (pre=col1, post=col2) with
%   dose_values = [pre_gv, post_gv], or [0, 1] as a placeholder.

    unique_ids  = unique(unit_ids_in_table, 'stable');
    n_units     = numel(unique_ids);
    unit_rec_ids = ctc.FeatureStore.UnitTable.RecordingID;

    if strcmp(method, 'metadata')
        detail = struct('fr_matrix', [], 'dose_values', [], 'unit_ids', {unique_ids});
        return
    end

    % Gather all unique GroupingVar values (sorted)
    if ismember(group_field, meta.Properties.VariableNames)
        all_gv = unique(meta.(group_field), 'sorted');
        if iscell(all_gv); all_gv = cell2mat(all_gv); end
        all_gv = all_gv(~isnan(all_gv));
    else
        all_gv = [0, 1];
    end
    n_doses = numel(all_gv);

    fr_matrix = nan(n_units, n_doses);

    for u = 1:n_units
        uid      = unique_ids(u);
        uid_rows = find(unit_ids_in_table == uid);
        if isempty(uid_rows); continue; end

        for ri = 1:numel(uid_rows)
            row    = uid_rows(ri);
            ud_idx = ud_order(row);
            if ud_idx <= 0; continue; end

            rec_id  = unit_rec_ids(row);
            rec_row = meta(meta.RecordingID == rec_id, :);
            if isempty(rec_row); continue; end

            if ismember(group_field, meta.Properties.VariableNames)
                gv = rec_row.(group_field)(1);
                if ~isnumeric(gv); gv = str2double(string(gv)); end
            else
                gv = NaN;
            end

            dose_col = find(all_gv == gv, 1);
            if isempty(dose_col); continue; end

            st = ud(ud_idx).SpikeTimes;
            dur = ud(ud_idx).RecordingDuration;
            fr_matrix(u, dose_col) = numel(st) / max(dur, eps);
        end
    end

    detail = struct('fr_matrix',  fr_matrix, ...
                    'dose_values', all_gv(:)', ...
                    'unit_ids',    {unique_ids});
end

%% ── Helper: per-unit full-curve dose-response ─────────────────────────────

function [flag_inc, flag_dec, strengths, pvals] = fullCurveResponsive( ...
        ud, ud_indices, unit_rows, unit_ids_in_table, ...
        meta, rec_ids_here, unit_rec_ids, unit_mask, ud_order, ...
        group_field, p)
% Three-stage responsive unit identification:
%   Stage 1 — Monotonicity: |Spearman rho(dose rank, FR)| >= threshold
%   Stage 2 — Effect size:  fold change (baseline to max dose) >= threshold
%   Stage 3 — Bootstrap:    permutation test on baseline vs max dose
%
% Vectorized where possible. Stage 3 only runs on units passing stages 1+2.

    direction = lower(string(p.Direction));
    n_rows    = numel(unit_rows);
    flag_inc  = false(1, n_rows);
    flag_dec  = false(1, n_rows);
    strengths = zeros(1, n_rows);
    pvals     = nan(1, n_rows);

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
        local_mask  = row_unit_ids == unique_uids(u);
        uid_frs     = all_frs(local_mask);
        uid_rec_ids = all_rec_ids(local_mask);
        for ri = 1:n_valid_recs
            hit = uid_rec_ids == valid_rec_str(ri);
            if any(hit)
                fr_matrix(u, ri) = uid_frs(find(hit, 1));
            end
        end
    end

    % ══════════════════════════════════════════════════════════════════════════
    % STAGE 1 — Monotonicity: Spearman rho(dose rank, FR) >= threshold
    % ══════════════════════════════════════════════════════════════════════════
    rhos = nan(n_unique, 1);

    full_mask = all(~isnan(fr_matrix), 2);
    if any(full_mask)
        fr_full = fr_matrix(full_mask, :)';
        r = corr(gv_valid, fr_full, 'Type', 'Spearman');
        rhos(full_mask) = r';
    end

    partial_idx = find(~full_mask & sum(~isnan(fr_matrix), 2) >= 2);
    for u = partial_idx'
        valid_pts = ~isnan(fr_matrix(u, :));
        try
            rhos(u) = corr(gv_valid(valid_pts), fr_matrix(u, valid_pts)', 'Type', 'Spearman');
        catch
        end
    end

    switch direction
        case "increase"
            pass_s1 = rhos >= p.MonotonicityThreshold;
        case "decrease"
            pass_s1 = rhos <= -p.MonotonicityThreshold;
        case "both"
            pass_s1 = abs(rhos) >= p.MonotonicityThreshold;
    end
    pass_s1(isnan(rhos)) = false;
    n_s1 = sum(pass_s1);

    % ══════════════════════════════════════════════════════════════════════════
    % STAGE 2 — Effect size: fold change (baseline to max dose) >= threshold
    % ══════════════════════════════════════════════════════════════════════════
    [~, baseline_col] = min(gv_valid);
    [~, maxdose_col]  = max(gv_valid);

    baseline_fr = fr_matrix(:, baseline_col);
    maxdose_fr  = fr_matrix(:, maxdose_col);

    % Fold change: ratio of higher/lower FR. Guard against zero FR.
    fr_floor    = 0.01;  % Hz — floor to avoid division by zero / infinite fold change
    base_safe   = max(baseline_fr, fr_floor);
    max_safe    = max(maxdose_fr,  fr_floor);

    switch direction
        case "increase"
            fold_change = max_safe ./ base_safe;
        case "decrease"
            fold_change = base_safe ./ max_safe;
        case "both"
            fold_change = max(max_safe ./ base_safe, base_safe ./ max_safe);
    end

    pass_s2 = pass_s1 & fold_change >= p.MinFoldChange & ...
              ~isnan(baseline_fr) & ~isnan(maxdose_fr);
    n_s2 = sum(pass_s2);

    % ══════════════════════════════════════════════════════════════════════════
    % STAGE 3 — Bootstrap: permutation test on baseline vs max dose
    % ══════════════════════════════════════════════════════════════════════════
    s3_idx = find(pass_s2);  % indices into unique_uids
    pass_s3 = false(n_unique, 1);
    boot_pvals = nan(n_unique, 1);

    if ~isempty(s3_idx)
        % Identify baseline and max-dose recording IDs
        baseline_rec = valid_rec_ids(baseline_col);
        maxdose_rec  = valid_rec_ids(maxdose_col);

        % Build matched pre/post UnitData indices for stage-3 units
        pre_ud_idx  = [];
        post_ud_idx = [];
        s3_kept     = [];  % which s3_idx entries have matched pre+post

        for si = 1:numel(s3_idx)
            u = s3_idx(si);
            uid = unique_uids(u);
            local_mask = row_unit_ids == uid;
            local_ud_raw = ud_order(unit_rows(local_mask));
            valid_ud = local_ud_raw > 0;
            local_ud = local_ud_raw(valid_ud);
            if isempty(local_ud); continue; end

            uid_rids = string({ud(local_ud).RecordingID})';
            pre_hit  = find(uid_rids == string(baseline_rec), 1);
            post_hit = find(uid_rids == string(maxdose_rec),  1);
            if isempty(pre_hit) || isempty(post_hit); continue; end

            pre_ud_idx  = [pre_ud_idx,  local_ud(pre_hit)];   %#ok<AGROW>
            post_ud_idx = [post_ud_idx, local_ud(post_hit)];  %#ok<AGROW>
            s3_kept     = [s3_kept, u];                        %#ok<AGROW>
        end

        if ~isempty(pre_ud_idx)
            % Build binned spike matrices
            pre_dur  = ud(pre_ud_idx(1)).RecordingDuration;
            post_dur = ud(post_ud_idx(1)).RecordingDuration;
            pre_win  = [0, pre_dur];
            post_win = [0, post_dur];

            n_bins   = floor(min(diff(pre_win), diff(post_win)) / p.BinSize);
            pre_win  = [0, n_bins * p.BinSize];
            post_win = [0, n_bins * p.BinSize];

            pre_mat  = CellTypeClassifier.binnedSpikeMatrix(ud(pre_ud_idx),  pre_win,  p.BinSize);
            post_mat = CellTypeClassifier.binnedSpikeMatrix(ud(post_ud_idx), post_win, p.BinSize);

            response = bootstrapFiringRateResponse(pre_mat, post_mat, ...
                p.ConfirmationNIter, p.ConfirmationAlpha);

            % Convert index vectors to logical masks
            n_tested = numel(s3_kept);
            is_inc = false(1, n_tested);
            is_dec = false(1, n_tested);
            is_inc(response.increase) = true;
            is_dec(response.decrease) = true;

            % Map bootstrap results back to unique-unit indices
            for ki = 1:n_tested
                u = s3_kept(ki);
                if isfield(response, 'p_values') && numel(response.p_values) >= ki
                    boot_pvals(u) = response.p_values(ki);
                end

                switch direction
                    case "increase"
                        pass_s3(u) = is_inc(ki);
                    case "decrease"
                        pass_s3(u) = is_dec(ki);
                    case "both"
                        pass_s3(u) = is_inc(ki) || is_dec(ki);
                end
            end
        end
    end
    n_s3 = sum(pass_s3);

    fprintf('  Full-curve pipeline: %d units -> S1 (monotonicity): %d -> S2 (fold change): %d -> S3 (bootstrap): %d\n', ...
        n_unique, n_s1, n_s2, n_s3);

    % ── Map results back to unit_rows ─────────────────────────────────────────
    for u = 1:n_unique
        local_idx            = find(row_unit_ids == unique_uids(u));
        strengths(local_idx) = abs(rhos(u));
        pvals(local_idx)     = boot_pvals(u);

        if pass_s3(u)
            if direction == "both"
                if rhos(u) > 0
                    flag_inc(local_idx) = true;
                else
                    flag_dec(local_idx) = true;
                end
            elseif direction == "increase"
                flag_inc(local_idx) = true;
            else
                flag_dec(local_idx) = true;
            end
        end
    end
end
