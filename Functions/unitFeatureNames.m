function names = unitFeatureNames(group)
% UNITFEATURENAMES  Column names of a single-unit feature group.
%
%   names = unitFeatureNames("act")   % activity features (computeActivityFeatures)
%   names = unitFeatureNames("reg")   % regularity features (computeActivityFeatures)
%   names = unitFeatureNames("bc")    % Bombcell quality metrics
%   names = unitFeatureNames("burst") % participation in network bursts (computeNetworkFeatures)
%
% Single source of the names used by computeActivityFeatures, the FeatureStore selectors
% and FeatureCatalog, so a feature added to the computation is also selected by its group.

arguments
    group (1,1) string {mustBeMember(group, ["act","reg","bc","burst"])}
end

switch group
    case "act"
        names = ["FiringRate","MeanInterSpikeInterval","VarianceInterSpikeInterval","CVInterSpikeInterval", ...
            "CV2InterSpikeInterval","LocalVariation","RevisedLocalVariation","FanoFactor","PartialAutocorrelation", ...
            "ISIBimodality"];
    case "reg"
        names = ["RegularityFrequency","RegularityMagnitude","RegularityFit"];
    case "burst"
        names = "BurstFractionSpikes";      % written by the network-burst step
    case "bc"
        names = ["BC_nPeaks","BC_nTroughs","BC_waveformDuration", ...
            "BC_mainPeakToTroughRatio","BC_peak1ToPeak2Ratio", ...
            "BC_scndPeakToTroughRatio","BC_troughToPeak2Ratio", ...
            "BC_waveformBaselineFlatness","BC_spatialDecaySlope", ...
            "BC_fractionRPVs","BC_presenceRatio","BC_percSpikesMissing"];
end
end
