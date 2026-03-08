function detectBursts(eia)
% DETECTBURSTS  Label each time bin as burst (2) or non-burst (1) state.
%
% A burst is defined as a time bin where:
%   - The I/E ratio is below the threshold (excitatory-dominant activity), AND
%   - The excitatory firing rate exceeds the noise floor
% The resulting binary vector is median-filtered to smooth rapid state transitions.
%
% Requires eia.Activity (run computeActivity first).
% Sets: eia.BurstState

arguments
    eia EIAnalyzer
end

assert(~isempty(eia.Activity), 'Run computeActivity() before detectBursts()');

p          = eia.Parameters.BurstDetection;
burst_state = cell(1, length(eia.Activity));

for c = 1:length(eia.Activity)
    act = eia.Activity(c);

    % Noise floor: second smallest unique value (avoids true zero from silent units)
    unique_E    = unique(act.E);
    noise_floor = unique_E(min(2, length(unique_E)));

    z = double(act.ratio < p.Threshold & act.E > noise_floor) + 1;
    z = medfilt1(z, p.SmoothWindow);
    burst_state{c} = z;
end

eia.BurstState = burst_state;
end
