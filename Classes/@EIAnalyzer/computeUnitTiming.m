function computeUnitTiming(eia)
% COMPUTEUNITTIMING  Every unit's firing relative to the network bursts of its recording
% (burstUnitTiming): latency of the first spike to onset, participation, rate before onset /
% in / after the burst relative to its mean rate. Labels (Parameters.Labels) are attached but
% not used, so E and I can be compared afterwards - also within firing-rate strata, as in
% summarizeBurstEI.
%
% Sets eia.UnitTiming (table): RecordingID, Culture, UnitRow (FeatureStore row), UnitID, Label,
% Unit, NSpikes, Rate, Latency, Participation, PreRatio, InRatio, PostRatio.
% Requires computeNetworkBursts.

arguments
    eia EIAnalyzer
end
assert(~isempty(eia.NetworkBursts), 'EIAnalyzer:noBursts', 'Run computeNetworkBursts() first.');
ud = eia.Classifier.UnitDataArray; map = eia.getCultureUnitMapping();
fs = eia.Classifier.FeatureStore; p = eia.Parameters.UnitTiming;
NB = eia.NetworkBursts; parts = cell(1, numel(NB));
for c = 1:numel(NB)
    rows = NB(c).UnitRows;
    if isempty(rows) || isempty(fieldnames(NB(c).Bursts)), continue; end
    u = ud(map.ud_order(rows));
    n_sp = arrayfun(@(x) numel(x.SpikeTimes), u);
    st = vertcat(u.SpikeTimes); su = repelem((1:numel(u))', n_sp(:));
    T = burstUnitTiming(st(:), su, NB(c).Bursts, NB(c).Duration, PreWindow=p.PreWindow, ...
        PostWindow=p.PostWindow, NUnits=numel(u));
    lab = map.labels(rows); lab(~ismember(lab, [1 2])) = NaN;
    T = [table(repmat(NB(c).RecordingID, numel(rows), 1), repmat(NB(c).Culture, numel(rows), 1), rows(:), ...
        string(fs.UnitTable.UnitID(rows)), lab(:), 'VariableNames', {'RecordingID', 'Culture', 'UnitRow', 'UnitID', 'Label'}), T]; %#ok<AGROW>
    parts{c} = T;
end
eia.UnitTiming = vertcat(parts{:});
end
