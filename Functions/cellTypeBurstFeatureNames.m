function [network_names, unit_names, metric_names] = cellTypeBurstFeatureNames()
% CELLTYPEBURSTFEATURENAMES  Column names written by computeCellTypeBurstEI.
%
%   [network_names, unit_names, metric_names] = cellTypeBurstFeatureNames()
%
% metric_names are the burstEIMetrics metrics (read from the function itself, so the list
% follows it); network_names = BurstEI_<metric>[_Excess] and BurstUnit<x>_E / _I / _IminusE /
% _IminusE_Excess; unit_names = BurstUnit<x>, x = the burstUnitTiming measures.

persistent cache
if isempty(cache)
    R = burstEIMetrics(zeros(0, 1), zeros(0, 1), zeros(1, 0), ...
        struct('Onset', [], 'Offset', [], 'T_peak', []), 1);
    m = string(R.Metrics.Metric)';
    x = ["Latency", "Participation", "PreRatio", "InRatio", "PostRatio"];
    ei = reshape(["BurstEI_" + m; "BurstEI_" + m + "_Excess"], 1, []);
    bu = reshape(["BurstUnit" + x + "_E"; "BurstUnit" + x + "_I"; ...
                  "BurstUnit" + x + "_IminusE"; "BurstUnit" + x + "_IminusE_Excess"], 1, []);
    cache = struct('net', [ei, bu], 'unit', "BurstUnit" + x, 'metric', m);
end
network_names = cache.net;
unit_names    = cache.unit;
metric_names  = cache.metric;
end
