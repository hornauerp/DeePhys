function [T, G] = summarizeBurstEI(eia, opts)
% SUMMARIZEBURSTEI  One row per recording: metadata, burst counts, every burstEIMetrics metric
% with its rate-matched null, and (when computed) unit timing and rhythm summaries.
%
%   T = eia.summarizeBurstEI()
%   [T, G] = eia.summarizeBurstEI(GroupBy=["DIV", "Construct", "Concentration"])
%
% Columns per metric m (see burstEIMetrics): m (value), m_Null (null median), m_PLower / m_PUpper
% (share of null draws <= / >= the value; e.g. IEContrastIn_PLower < 0.05 = I/E lower during
% bursts than for 95 % of random I sets with the same firing rates).
% Unit timing (computeUnitTiming): median Latency / Participation / PreRatio of I and E units,
% the rank-sum p of the latency difference, and LatencyDiffStrata - the median over firing-rate
% quintiles (of the recording's labelled units) of the I - E latency difference.
% Rhythm (computeBurstRhythm): Locked reliability and frequency of All / E / I; shared
% coherency of random halves of all units (bands) and its peak P.
% G (with GroupBy): per group the number of recordings, medians of the main metrics and the
% share of recordings with IEContrastIn_PLower < 0.05.

arguments
    eia EIAnalyzer
    opts.GroupBy string = string.empty
end
assert(~isempty(eia.BurstEI), 'EIAnalyzer:noBurstEI', 'Run computeBurstEI() first.');
meta = eia.Classifier.FeatureStore.MetadataTable;
B = eia.BurstEI; n = numel(B);
rows = cell(n, 1);
for c = 1:n
    r = table(B(c).RecordingID, B(c).Culture, 'VariableNames', {'RecordingID', 'Culture'});
    mrow = meta(string(meta.RecordingID) == B(c).RecordingID, :);
    for g = opts.GroupBy
        v = NaN; if ~isempty(mrow) && ismember(g, mrow.Properties.VariableNames), v = mrow.(g)(1); end
        if iscell(v), v = string(v); end
        r.(g) = v;
    end
    R = B(c).Result;
    if isempty(fieldnames(R)), rows{c} = []; continue; end
    r.Valid = R.Valid; r.NBursts = R.NBursts; r.NI = R.NI; r.NE = R.NE; r.NInnerBursts = R.NInnerBursts;
    M = R.Metrics;
    for k = 1:height(M)
        m = M.Metric(k);
        r.(m) = M.Value(k); r.(m + "_Null") = M.NullMedian(k);
        r.(m + "_PLower") = M.PLower(k); r.(m + "_PUpper") = M.PUpper(k);
    end
    r = addTiming(r, eia.UnitTiming, B(c).RecordingID);
    r = addRhythm(r, eia.BurstRhythm, B(c).RecordingID);
    rows{c} = r;
end
rows = rows(~cellfun(@isempty, rows));
T = vertcat(rows{:});
G = table();
if ~isempty(opts.GroupBy) && ~isempty(T)
    keyv = ["IEContrastIn", "IEContrastPre", "IEContrastPost", "IEMaxTime", "IEMinTime", "ITerminal", "AntiPhase"];
    if ismember("LatencyI", T.Properties.VariableNames), keyv = [keyv, "LatencyI", "LatencyE", "LatencyDiffStrata"]; end
    V = T(T.Valid, :);
    V.SigIEContrastIn = V.IEContrastIn_PLower < 0.05;
    G = groupsummary(V, cellstr(opts.GroupBy), {'median', 'mean'}, cellstr([keyv, "SigIEContrastIn"]));
    keep = startsWith(G.Properties.VariableNames, 'median_') | ismember(G.Properties.VariableNames, [cellstr(opts.GroupBy), {'GroupCount', 'mean_SigIEContrastIn'}]);
    G = G(:, keep);
end
end

function r = addTiming(r, U, rec)
if isempty(U), return; end
u = U(U.RecordingID == rec, :);
I = u(u.Label == 2, :); E = u(u.Label == 1, :);
r.LatencyI = median(I.Latency, 'omitnan'); r.LatencyE = median(E.Latency, 'omitnan');
r.ParticipationI = median(I.Participation, 'omitnan'); r.ParticipationE = median(E.Participation, 'omitnan');
r.PreRatioI = median(I.PreRatio, 'omitnan'); r.PreRatioE = median(E.PreRatio, 'omitnan');
r.LatencyP = NaN; r.LatencyDiffStrata = NaN;
if sum(isfinite(I.Latency)) >= 3 && sum(isfinite(E.Latency)) >= 3
    r.LatencyP = ranksum(I.Latency, E.Latency);
    lab = u(ismember(u.Label, [1 2]), :);
    q = prctile(lab.Rate, 20:20:80); q = q(:)';
    s = discretize(lab.Rate, unique([-inf, q, inf]));
    d = nan(1, max(s));
    for k = 1:max(s)
        a = lab.Latency(s == k & lab.Label == 2); b = lab.Latency(s == k & lab.Label == 1);
        if sum(isfinite(a)) >= 2 && sum(isfinite(b)) >= 2, d(k) = median(a, 'omitnan') - median(b, 'omitnan'); end
    end
    r.LatencyDiffStrata = median(d, 'omitnan');
end
end

function r = addRhythm(r, BR, rec)
if isempty(BR), return; end
for f = ["LockedReliabilityAll", "LockedPAll", "LockedFreqAll", "LockedReliabilityE", "LockedPE", "LockedFreqE", ...
        "LockedReliabilityI", "LockedPI", "LockedFreqI", "SharedCoh5_15", "SharedCoh15_30", "SharedCoh30_60", ...
        "SharedPeakFreq", "SharedPeakP", "SharedLagEI"]
    r.(f) = NaN;
end
k = find([BR.RecordingID] == rec, 1);
if isempty(k) || isempty(fieldnames(BR(k).Result)), return; end
L = BR(k).Result.Locked; S = BR(k).Result.Shared;
for g = ["All", "E", "I"]
    j = find(L.Group == g, 1);
    if isempty(j), continue; end
    r.("LockedReliability" + g) = L.Reliability(j); r.("LockedP" + g) = L.P(j); r.("LockedFreq" + g) = L.Freq(j);
end
j = find(S.Pair == "All halves", 1);
if ~isempty(j)
    r.SharedCoh5_15 = S.Coh5_15(j); r.SharedCoh15_30 = S.Coh15_30(j); r.SharedCoh30_60 = S.Coh30_60(j);
    r.SharedPeakFreq = S.PeakFreq(j); r.SharedPeakP = S.P(j);
end
j = find(S.Pair == "E-I", 1);
if ~isempty(j), r.SharedLagEI = S.Lag(j); end
end
