function acg = getACG(unit, bin_size, lag)
% GETACG  Recompute autocorrelogram from spike times at arbitrary parameters.
%
% Always recomputes from spike times — does NOT read unit.ACG or unit.FullACG.
% Use unit.FullACG directly when you want the pre-computed high-resolution ACG,
% and unit.ACG for the stored standard-resolution ACG.
%
% This function is intended for harmonization: when DeePhys FullACG parameters
% (0.5 ms / 100 ms) must match a different external dataset configuration.
%
% INPUTS:
%   unit     - Unit object
%   bin_size - bin width in seconds
%   lag      - one-sided lag in seconds
% OUTPUT:
%   acg - (n_bins x 1) double, raw spike counts (central peak intact)

arguments
    unit     Unit
    bin_size (1,1) double
    lag      (1,1) double
end

n_bins_expected = round(2 * lag / bin_size) + 1;

if isempty(unit.SpikeTimes)
    acg = zeros(n_bins_expected, 1);
    return
end

[acg_3d, ~] = CCG(unit.SpikeTimes, ones(size(unit.SpikeTimes)), ...
    'binSize', bin_size, ...
    'duration', 2 * lag, ...
    'Fs', 1 / unit.SamplingRate);
acg = double(acg_3d(:, 1, 1));
end
