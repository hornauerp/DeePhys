function [spike_times, spike_units, truth] = simulateNetworkBursts(opts)
% SIMULATENETWORKBURSTS  Synthetic population spike trains with known network bursts.
%
%   [spike_times, spike_units, truth] = simulateNetworkBursts()
%   [...] = simulateNetworkBursts(Mode="independent", NTonic=10, SubBursts=3)
%
% Ground truth for burst detectors and the E/I analysis. Every unit fires Poisson background
% spikes; tonic units fire regularly (gamma renewal, CV ~0.3) and ignore the bursts. Bursting
% units fire at InBurstRate inside bursts they take part in (probability Participation per
% burst, lognormal gain per unit).
%   Mode "network"     - all bursting units share the burst times (network bursts)
%   Mode "independent" - every unit bursts on its own schedule with the same statistics:
%                        single-unit bursts but no network synchrony (negative control)
%   Mode "none"        - background and tonic spikes only (negative control)
% SubBursts > 1 splits each burst into evenly spaced sub-events (reverberations) covering
% SubDutyCycle of the burst, with the same expected spike count.
%
% OUTPUTS:
%   spike_times - (S x 1) sorted spike times (s)
%   spike_units - (S x 1) unit 1..NUnits
%   truth       - struct: T_start, T_end (network bursts, empty unless "network"),
%                 Sub_start, Sub_end (sub-events), Tonic (1 x NUnits logical), Mode
arguments
    opts.NUnits (1,1) double = 60
    opts.Duration (1,1) double = 300
    opts.BackgroundRate (1,1) double = 0.5      % Hz, median over units (lognormal, sigma 0.5)
    opts.NTonic (1,1) double = 6
    opts.TonicRate (1,2) double = [5 20]        % Hz, uniform range
    opts.BurstRate (1,1) double = 0.1           % bursts per second
    opts.BurstDuration (1,1) double = 0.3       % s, mean (gamma, CV 0.3)
    opts.InBurstRate (1,1) double = 30          % Hz per participating unit
    opts.Participation (1,1) double = 0.8
    opts.SubBursts (1,1) double = 1
    opts.SubDutyCycle (1,1) double = 0.5
    opts.Mode (1,1) string {mustBeMember(opts.Mode, ["network", "independent", "none"])} = "network"
    opts.Seed (1,1) double = 1
end
prev = rng; rng(opts.Seed, 'twister'); restore = onCleanup(@() rng(prev));   % global stream: poissrnd / gamrnd use it
U = opts.NUnits; T = opts.Duration;
tonic = false(1, U); tonic(U - opts.NTonic + 1:U) = true;
st = cell(U, 1);
% background
bg = opts.BackgroundRate * exp(0.5 * randn(U, 1));
for u = 1:U
    st{u} = T * rand(poissrnd(bg(u) * T), 1);
end
% tonic units: gamma renewal (shape 10)
for u = find(tonic)
    r = opts.TonicRate(1) + diff(opts.TonicRate) * rand;
    isi = gamrnd(10, 1 / (10 * r), ceil(1.5 * r * T) + 10, 1);
    t = cumsum(isi) + rand / r;
    st{u} = [st{u}; t(t < T)];
end
% bursts
truth = struct('T_start', zeros(1, 0), 'T_end', zeros(1, 0), 'Sub_start', zeros(1, 0), ...
    'Sub_end', zeros(1, 0), 'Tonic', tonic, 'Mode', opts.Mode);
bursting = find(~tonic);
gain = exp(0.4 * randn(U, 1));
if opts.Mode == "network"
    [b0, b1] = burstSchedule(opts);
    [s0, s1] = subEvents(b0, b1, opts);
    truth.T_start = b0; truth.T_end = b1; truth.Sub_start = s0; truth.Sub_end = s1;
    for u = bursting
        st{u} = [st{u}; burstSpikes(s0, s1, b0, b1, gain(u), opts)];
    end
elseif opts.Mode == "independent"
    for u = bursting
        [b0, b1] = burstSchedule(opts);
        [s0, s1] = subEvents(b0, b1, opts);
        st{u} = [st{u}; burstSpikes(s0, s1, b0, b1, gain(u), opts)];
    end
end
n = cellfun(@numel, st);
spike_times = vertcat(st{:});
spike_units = repelem((1:U)', n);
[spike_times, o] = sort(spike_times);
spike_units = spike_units(o);
end

function [b0, b1] = burstSchedule(opts)
% burst onsets: exponential intervals plus a refractory gap of two burst durations
b0 = zeros(1, 0); b1 = b0;
if opts.BurstRate <= 0, return; end
t = exprnd(1 / opts.BurstRate);
while t < opts.Duration - 1
    d = gamrnd(1 / 0.3^2, opts.BurstDuration * 0.3^2);
    b0(end+1) = t; b1(end+1) = min(t + d, opts.Duration); %#ok<AGROW>
    t = t + d + 2 * opts.BurstDuration + exprnd(1 / opts.BurstRate);
end
end

function [s0, s1] = subEvents(b0, b1, opts)
k = max(1, round(opts.SubBursts));
if k == 1, s0 = b0; s1 = b1; return; end
s0 = zeros(1, 0); s1 = s0;
for i = 1:numel(b0)
    d = b1(i) - b0(i); w = d * opts.SubDutyCycle / k; gap = (d - k * w) / max(k - 1, 1);
    on = b0(i) + (0:k-1) * (w + gap);
    s0 = [s0, on]; s1 = [s1, on + w]; %#ok<AGROW>
end
end

function t = burstSpikes(s0, s1, b0, b1, g, opts)
% Poisson spikes in the sub-events of the bursts this unit takes part in; the rate is scaled
% so that a burst holds InBurstRate x burst duration spikes on average
k = max(1, round(opts.SubBursts));
duty = 1; if k > 1, duty = opts.SubDutyCycle; end
part = rand(numel(b0), 1) < opts.Participation;
t = zeros(0, 1);
for i = find(part)'
    in = s0 >= b0(i) & s1 <= b1(i) + 1e-12;
    for j = find(in)
        n = poissrnd(g * opts.InBurstRate / duty * (s1(j) - s0(j)));
        t = [t; s0(j) + (s1(j) - s0(j)) * rand(n, 1)]; %#ok<AGROW>
    end
end
end
