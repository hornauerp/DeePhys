function identifyResponsiveUnits(ctc, metadata_filter)
<<<<<<< HEAD
% IDENTIFYRESPONSIVEUNITS  Bootstrap test to identify units with significant firing rate changes.
%
% Groups units by culture (identity keys in MetadataTable), then for each
% culture compares pre/post spike rates with a bootstrap permutation test.
% Units showing a significant change in the configured direction are flagged
% as responsive candidates.
=======
% IDENTIFYRESPONSIVEUNITS  Identify units with significant firing rate changes.
%
% Supports two ground truth methods (set via Parameters.Bootstrap.GroundTruthMethod):
%
%   "full_curve" (default) — For each culture, computes mean firing rate per
%       recording (concentration level) and tests for monotonic increase using
%       Spearman rank correlation with a permutation test. Cultures with fewer
%       than MinRecordings recordings fall back to the two-window method.
%       Strength = absolute Spearman rho (0-1).
%
%   "two_window" — Legacy method comparing pre-treatment vs post-treatment
%       firing rates using a bootstrap permutation test.
%       Strength = |empirical difference| / null std (effect size).
%
% Units that significantly increase are flagged as inhibitory candidates
% (ResponsiveUnitIdx = true). For each unit, a continuous strength score is
% stored in ctc.ResponsiveStrength — higher values indicate more credible
% responsive identification.
%
% When MinResponsiveStrength > 0, responsive units below this strength
% threshold are demoted to non-responsive. This filters out units that
% barely pass the significance test but have weak dose-response trends.
>>>>>>> claude/dreamy-varahamihira
%
% INPUTS:
%   ctc             - CellTypeClassifier (new API with FeatureStore + UnitDataArray)
%   metadata_filter - (optional) {field, value} cell pair to restrict which cultures
<<<<<<< HEAD
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
=======
%                     are tested (e.g. {'AAV', 128}). Default: all cultures.
%
% Sets: ctc.ResponsiveUnitIdx, ctc.ResponsiveStrength, ctc.UnitList
>>>>>>> claude/dreamy-varahamihira

arguments
    ctc             CellTypeClassifier
    metadata_filter cell = {}
end

<<<<<<< HEAD
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
=======
p = ctc.Parameters.Bootstrap;
rg = ctc.RecordingGroup;
method = p.GroundTruthMethod;

responsive_idx = [];
strength_all   = [];
unit_list = Unit.empty;
>>>>>>> claude/dreamy-varahamihira

% Build culture ID per unit (from MetadataTable via RecordingID join)
meta            = fs.MetadataTable;
unit_rec_ids    = fs.UnitTable.RecordingID;
culture_ids     = FeatureStore.getCultureIDsForUnits( ...
    unit_rec_ids, meta, ctc.Parameters.CultureKeys);
unique_cultures = unique(culture_ids, 'stable');

<<<<<<< HEAD
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
=======
    % Apply optional metadata filter — only count increases from matching cultures
    if ~isempty(metadata_filter) && isfield(culture.Metadata, metadata_filter{1})
        culture_matches = any([culture.Metadata.(metadata_filter{1})] == metadata_filter{2});
    else
        culture_matches = true;
    end

    % Select ground truth method
    if method == "full_curve" && culture.N >= p.MinRecordings
        response = culture.bootstrapDoseResponse(p.NIter, p.Alpha);
    else
        if method == "full_curve" && culture.N < p.MinRecordings
            fprintf('Culture %d has only %d recordings (<%d) — falling back to two_window\n', ...
                c, culture.N, p.MinRecordings);
        end
        response = culture.bootstrapResponse(p.PreCutout, p.PostCutout, p.BinSize, p.NIter, p.Alpha);
    end

    n_prev = length(unit_list);
    if culture_matches
        responsive_idx = [responsive_idx, response.increase + n_prev]; %#ok<AGROW>
    end
    strength_all = [strength_all, response.strength]; %#ok<AGROW>
    unit_list = [unit_list, culture.Units]; %#ok<AGROW>
end

ctc.UnitList = unit_list;
ctc.ResponsiveUnitIdx = false(1, length(unit_list));
ctc.ResponsiveUnitIdx(responsive_idx) = true;
ctc.ResponsiveStrength = strength_all;

% ── Optional strength-based filtering ────────────────────────────────────────
if isfield(p, 'MinResponsiveStrength') && p.MinResponsiveStrength > 0
    weak = ctc.ResponsiveUnitIdx & (strength_all < p.MinResponsiveStrength);
    n_demoted = sum(weak);
    ctc.ResponsiveUnitIdx(weak) = false;
    if n_demoted > 0
        fprintf('Demoted %d weak responders (strength < %.2f)\n', ...
            n_demoted, p.MinResponsiveStrength);
    end
end

fprintf('Inhibitory candidates: %i / %i units (%.1f%%)\n', ...
    sum(ctc.ResponsiveUnitIdx), length(ctc.UnitList), ...
    100 * mean(ctc.ResponsiveUnitIdx));
>>>>>>> claude/dreamy-varahamihira
end
