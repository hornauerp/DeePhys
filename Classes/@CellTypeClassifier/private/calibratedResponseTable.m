function T = calibratedResponseTable(ud_all, meta, culture_keys, group_field, baseline_gv, cultures_sel, p)
% CALIBRATEDRESPONSETABLE  Per-unit drug response (baseline vs highest doses), drift-corrected.
%
% For every unit present at baseline and at each of the p.HighDoses highest finite doses of
% its culture:
%   LFC  = log2 fold change of the firing rate (high doses pooled vs baseline; half-spike
%          pseudo-count)
%   LFCc = LFC minus the culture's median LFC (removes culture-wide drift / network effects)
%   p    = rank-sum test of p.BinSize-s spike counts (high doses vs baseline)
% Thresholds are set elsewhere from control cultures (control_calibrated method).
%
% INPUTS:
%   ud_all       - UnitData of all recordings (before discardSegmentUnits)
%   meta         - FeatureStore.MetadataTable (all recordings)
%   culture_keys - Parameters.CultureKeys
%   group_field  - Parameters.UMAP.GroupingVar (dose column)
%   baseline_gv  - Parameters.UMAP.GroupingValues (baseline dose)
%   cultures_sel - (string) culture IDs to evaluate
%   p            - Parameters.Bootstrap (BinSize, HighDoses)
% OUTPUT:
%   T - table: UnitID, Culture, BaselineRecordingID, FR0, FRhigh, LFC, LFCc, p
    rec_ids  = string(meta.RecordingID);
    rec_cult = FeatureStore.buildCultureIDs(meta, culture_keys);
    gv = meta.(group_field);
    if ~isnumeric(gv), gv = str2double(string(gv)); end
    ud_rec = string({ud_all.RecordingID});
    ud_uid = string({ud_all.UnitID});
    n_high = 2;
    if isfield(p, 'HighDoses'), n_high = p.HighDoses; end
    bin = p.BinSize;

    rows = cell(0, 8);
    for c = cultures_sel(:)'
        rmask = rec_cult == c & isfinite(gv);
        r_here = rec_ids(rmask); g_here = gv(rmask);
        base = r_here(g_here == baseline_gv);
        doses = sort(unique(g_here(g_here ~= baseline_gv & g_here > baseline_gv)));
        if isempty(base) || numel(doses) < n_high, continue, end
        base = base(1);
        high = r_here(ismember(g_here, doses(end-n_high+1:end)));
        bi = find(ud_rec == base);
        hi_sets = arrayfun(@(r) find(ud_rec == r), high, 'UniformOutput', false);
        lfc = nan(numel(bi), 1); pv = lfc; r0 = lfc; rh = lfc; keep = false(numel(bi), 1);
        for k = 1:numel(bi)
            u = ud_all(bi(k));
            cb = histcounts(u.SpikeTimes, 0:bin:u.RecordingDuration);
            ch = [];
            ok = true;
            for h = 1:numel(hi_sets)
                j = hi_sets{h}(ud_uid(hi_sets{h}) == string(u.UnitID));
                if isempty(j), ok = false; break, end
                v = ud_all(j(1));
                ch = [ch, histcounts(v.SpikeTimes, 0:bin:v.RecordingDuration)]; %#ok<AGROW>
            end
            if ~ok || isempty(cb) || isempty(ch), continue, end
            keep(k) = true;
            r0(k) = mean(cb) / bin; rh(k) = mean(ch) / bin;
            lfc(k) = log2((rh(k) + 0.5 / (numel(ch) * bin)) / (r0(k) + 0.5 / (numel(cb) * bin)));
            pv(k) = ranksum(ch, cb);
        end
        if ~any(keep), continue, end
        lfcc = lfc(keep) - median(lfc(keep));
        n = sum(keep);
        rows(end+1, :) = {ud_uid(bi(keep))', repmat(c, n, 1), repmat(base, n, 1), ...
            r0(keep), rh(keep), lfc(keep), lfcc, pv(keep)}; %#ok<AGROW>
    end
    if isempty(rows)
        T = table(strings(0, 1), strings(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
            'VariableNames', {'UnitID', 'Culture', 'BaselineRecordingID', 'FR0', 'FRhigh', 'LFC', 'LFCc', 'p'});
        return
    end
    T = table(vertcat(rows{:, 1}), vertcat(rows{:, 2}), vertcat(rows{:, 3}), vertcat(rows{:, 4}), ...
        vertcat(rows{:, 5}), vertcat(rows{:, 6}), vertcat(rows{:, 7}), vertcat(rows{:, 8}), ...
        'VariableNames', {'UnitID', 'Culture', 'BaselineRecordingID', 'FR0', 'FRhigh', 'LFC', 'LFCc', 'p'});
end
