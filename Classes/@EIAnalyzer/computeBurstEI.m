function computeBurstEI(eia)
% COMPUTEBURSTEI  E / I activity around the network bursts of every recording (burstEIMetrics).
%
% Labels from Parameters.Labels (graph / composition / confident / responders); unlabelled units
% (NaN) are left out of both populations but still took part in finding the bursts. Options:
% Parameters.BurstEI (PreWindow, PostWindow, Margin, MinBetween, BinSize, Window, NNull, NullStrata,
% MinUnits, MinBursts, Seed). Recordings with too few bursts or units get Result.Valid = false.
%
% Sets eia.BurstEI (1 x N_recordings): RecordingID, Culture, Labels (of NetworkBursts.UnitRows),
% Result. Requires computeNetworkBursts.

arguments
    eia EIAnalyzer
end
assert(~isempty(eia.NetworkBursts), 'EIAnalyzer:noBursts', 'Run computeNetworkBursts() first.');
ud = eia.Classifier.UnitDataArray; map = eia.getCultureUnitMapping();
labels_all = map.labels;
p = eia.Parameters.BurstEI;
opts = {'PreWindow', p.PreWindow, 'PostWindow', p.PostWindow, 'Margin', p.Margin, 'MinBetween', p.MinBetween, 'BinSize', p.BinSize, ...
    'Window', p.Window, 'NNull', p.NNull, 'NullStrata', p.NullStrata, 'MinUnits', p.MinUnits, ...
    'MinBursts', p.MinBursts, 'Seed', p.Seed};
NB = eia.NetworkBursts; n = numel(NB);
out = repmat(struct('RecordingID', "", 'Culture', "", 'Labels', [], 'Result', struct()), 1, n);
for c = 1:n
    out(c).RecordingID = NB(c).RecordingID; out(c).Culture = NB(c).Culture;
    rows = NB(c).UnitRows;
    if isempty(rows) || isempty(fieldnames(NB(c).Bursts)), continue; end
    lab = labels_all(rows); lab(~ismember(lab, [1 2])) = NaN;
    out(c).Labels = lab;
    [st, su] = spikesOf(ud(map.ud_order(rows)));
    out(c).Result = burstEIMetrics(st, su, lab, NB(c).Bursts, NB(c).Duration, opts{:});
end
eia.BurstEI = out;
end

function [st, su] = spikesOf(u)
n_sp = arrayfun(@(x) numel(x.SpikeTimes), u);
st = vertcat(u.SpikeTimes); st = st(:);
su = repelem((1:numel(u))', n_sp(:));
end
