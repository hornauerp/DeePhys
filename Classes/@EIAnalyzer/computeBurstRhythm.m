function computeBurstRhythm(eia)
% COMPUTEBURSTRHYTHM  Fast structure inside the network bursts of every recording (burstRhythm):
% split-half reliability and frequency of the burst-locked sub-peaks of E, I and all units, and
% shared fluctuations that are not locked to the burst (random halves of all units, and E vs I
% with the lag of I in the cycle), each against its null. Options: Parameters.BurstRhythm.
%
% Sets eia.BurstRhythm (1 x N_recordings): RecordingID, Culture, Result. Requires
% computeNetworkBursts.

arguments
    eia EIAnalyzer
end
assert(~isempty(eia.NetworkBursts), 'EIAnalyzer:noBursts', 'Run computeNetworkBursts() first.');
ud = eia.Classifier.UnitDataArray; map = eia.getCultureUnitMapping();
p = eia.Parameters.BurstRhythm;
NB = eia.NetworkBursts; n = numel(NB);
out = repmat(struct('RecordingID', "", 'Culture', "", 'Result', struct()), 1, n);
for c = 1:n
    out(c).RecordingID = NB(c).RecordingID; out(c).Culture = NB(c).Culture;
    rows = NB(c).UnitRows;
    if isempty(rows) || isempty(fieldnames(NB(c).Bursts)), continue; end
    u = ud(map.ud_order(rows));
    n_sp = arrayfun(@(x) numel(x.SpikeTimes), u);
    st = vertcat(u.SpikeTimes); su = repelem((1:numel(u))', n_sp(:));
    lab = map.labels(rows);
    groups = struct('All', true(1, numel(rows)), 'E', lab == 1, 'I', lab == 2);
    out(c).Result = burstRhythm(st(:), su, groups, NB(c).Bursts, BinSize=p.BinSize, SmoothSigma=p.SmoothSigma, ...
        LockedWindow=p.LockedWindow, SharedWindow=p.SharedWindow, NSplits=p.NSplits, NShift=p.NShift, ...
        MinBursts=p.MinBursts, Pairs={"All"; "E"; ["E", "I"]});
end
eia.BurstRhythm = out;
end
