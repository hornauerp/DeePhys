function detectBursts(eia)
% DETECTBURSTS  Label each time bin as burst (2) or non-burst (1) state.
%
% A burst is defined as a time bin where:
%   - The I/E ratio is below MaxIERatio (excitatory-dominant activity), AND
%   - The excitatory firing rate exceeds the noise floor
% The resulting binary vector is median-filtered to smooth rapid state transitions.
%
% Noise floor: NoiseFloorPercentile-th percentile of positive E values.
% This is robust to the continuous-valued E trace produced by movmean smoothing,
% unlike the "second unique value" heuristic which degenerates for float data.
%
% Clears BurstCutouts and NormalizedCutouts to prevent stale results after
% re-running with new parameters.
%
% Requires eia.Activity (run computeActivity first).
% Sets: eia.BurstState

arguments
    eia EIAnalyzer
end

assert(~isempty(eia.Activity), 'Run computeActivity() before detectBursts()');

% Invalidate downstream results
eia.BurstCutouts      = [];
eia.NormalizedCutouts = [];

p          = eia.Parameters.BurstDetection;
burst_state = cell(1, numel(eia.Activity));

for c = 1:numel(eia.Activity)
    act = eia.Activity(c);

    if isempty(act.E)
        burst_state{c} = [];
        continue
    end

    % Percentile-based noise floor over positive E values
    pos_E = act.E(act.E > 0);
    if isempty(pos_E)
        noise_floor = 0;
    else
        noise_floor = prctile(pos_E, p.NoiseFloorPercentile);
    end

    z = double(act.ratio < p.MaxIERatio & act.E > noise_floor) + 1;
    z = round(medfilt1(z, p.SmoothWindow));
    burst_state{c} = z;
end

eia.BurstState = burst_state;
end
