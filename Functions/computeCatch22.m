function catch22_table = computeCatch22(spike_times, duration, bin_size)
% COMPUTECATCH22  Compute catch22 features from a spike train.
%
%   catch22_table = computeCatch22(spike_times, duration)
%   catch22_table = computeCatch22(spike_times, duration, bin_size)
%
% INPUTS:
%   spike_times - (N x 1) double, spike times in seconds
%   duration    - (1 x 1) double, recording duration in seconds
%   bin_size    - (1 x 1) double, bin width in seconds (default 0.1)
%
% OUTPUTS:
%   catch22_table - (1 x 22) table of catch22 feature values

    arguments
        spike_times (:,1) double
        duration    (1,1) double
        bin_size    (1,1) double = 0.1
    end

    n_bins = round(duration / bin_size);
    binned = histcounts(spike_times, n_bins);
    [featureValues, featureNames] = catch22_all(binned');
    catch22_table = array2table(featureValues', 'VariableNames', featureNames);
end
