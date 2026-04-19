function detectBursts(eia)
% DETECTBURSTS  Label each time bin as burst (true) or non-burst (false).
%
% A burst is defined as a time bin where:
%   - The I/E ratio is below MaxIERatio (excitatory-dominant activity), AND
%   - The excitatory firing rate exceeds the noise floor
% The resulting logical vector is median-filtered to smooth rapid state transitions.
%
% Noise floor: NoiseFloorPercentile-th percentile of positive E values.
% This is robust to the continuous-valued E trace produced by movmean smoothing,
% unlike the "second unique value" heuristic which degenerates for float data.
% Note: the noise floor is computed from the smoothed E trace (Activity.E), so it
% implicitly depends on Activity.SmoothWindow.
%
% Clears BurstCutouts and NormalizedCutouts to prevent stale results after
% re-running with new parameters.
%
% Requires eia.Activity (run computeActivity first).
% Sets: eia.BurstState  (logical vector per culture: true = burst bin)

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
    fprintf('  detectBursts: culture %d/%d\n', c, numel(eia.Activity));
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

    z = act.ratio < p.MaxIERatio & act.E > noise_floor;
    z = logical(round(medfilt1(double(z), p.SmoothWindow)));
    burst_state{c} = z;
end

eia.BurstState = burst_state;
end
