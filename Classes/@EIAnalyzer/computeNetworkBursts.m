function computeNetworkBursts(eia)
% COMPUTENETWORKBURSTS  Network bursts of every recording with detectNetworkBursts.
%
% All units of a recording in the classifier take part (whatever their label), so the bursts do
% not depend on the E / I labels. Options of the detector: Parameters.NetworkBursts.Options
% (struct of name-value pairs, e.g. struct('MinQuietGap', 0.3)). Duration = the units'
% RecordingDuration (else the last spike).
%
% Sets eia.NetworkBursts (1 x N_recordings): RecordingID, Culture, UnitRows (FeatureStore rows,
% unit k of the recording = UnitRows(k)), Duration, Bursts, Info. Clears BurstEI, UnitTiming
% and BurstRhythm.

arguments
    eia EIAnalyzer
end
map = eia.getCultureUnitMapping();
ud = eia.Classifier.UnitDataArray;
opts = namedargs(eia.Parameters.NetworkBursts.Options);
n = numel(map.unique_cultures);
nb = repmat(struct('RecordingID', "", 'Culture', "", 'UnitRows', [], 'Duration', NaN, ...
    'Bursts', struct(), 'Info', struct()), 1, n);
for c = 1:n
    rows = find(map.culture_ids == map.unique_cultures(c) & map.ud_order > 0);
    nb(c).RecordingID = map.unique_cultures(c); nb(c).Culture = map.culture_names(c); nb(c).UnitRows = rows(:)';
    if isempty(rows), continue; end
    u = ud(map.ud_order(rows));
    [st, su] = spikesOf(u);
    dur = max([u.RecordingDuration]);
    if isempty(dur) || ~isfinite(dur), dur = max(st) + eps; end
    nb(c).Duration = dur;
    fprintf('  computeNetworkBursts: recording %d/%d (%s), %d units\n', c, n, map.unique_cultures(c), numel(rows));
    [nb(c).Bursts, nb(c).Info] = detectNetworkBursts(st, su, dur, opts{:});
end
eia.NetworkBursts = nb;
eia.BurstEI = []; eia.UnitTiming = []; eia.BurstRhythm = [];
end

function [st, su] = spikesOf(u)
n_sp = arrayfun(@(x) numel(x.SpikeTimes), u);
st = vertcat(u.SpikeTimes); st = st(:);
su = repelem((1:numel(u))', n_sp(:));
end

function c = namedargs(s)
f = fieldnames(s); c = cell(1, 2 * numel(f));
c(1:2:end) = f; c(2:2:end) = struct2cell(s);
end
