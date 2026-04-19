classdef EIAnalyzer < handle
% EIANALYZER  Excitatory/inhibitory network activity analysis.
%
% Takes a trained CellTypeClassifier and computes per-culture E/I activity
% traces, burst detection, aligned burst cutout extraction, and unit-population
% correlations.
%
% USAGE:
%   eia = EIAnalyzer(ctc, params);
%   eia.computeActivity();        % E, I, total firing rate traces per culture
%   eia.detectBursts();           % binary burst-state labelling from I/E ratio
%   eia.extractBurstCutouts();    % outlier-rejected + clustered burst cutouts
%   eia.normalizeBurstCutouts();  % cross-culture normalised burst profiles
%   eia.computeCorrelations();    % unit <-> population E/I correlations
%
% Plotting:
%   eia.PlotNetworkActivity(c)    % firing rate coloured by burst state
%   eia.PlotBurstCutouts(c)       % mean E/I burst shapes and gradients
%   eia.PlotEITraces(c)           % stacked E/I + total firing rate
%   eia.PlotBurstRaster(c)        % unit raster aligned to burst peaks (E=blue, I=red)
%
% BurstCutouts structure (per culture c):
%   .all           — all bursts after initial xcorr alignment
%   .accepted      — bursts surviving outlier rejection
%   .rejected      — discarded outlier bursts
%   .clusters      — (1 x K) struct array; each cluster xcorr-aligned to its own mean
%   .cluster_idx   — assignment of each accepted burst to a cluster (1 x N_accepted)
%   .silhouette    — mean silhouette score for the chosen k (quality metric)
%   .accepted_locs — (1 x N_accepted) burst peak bin indices into Activity.total
%                    (used for spike-time alignment in PlotBurstRaster)

    properties
        Classifier      % CellTypeClassifier — provides UnitLabels, FeatureStore, UnitDataArray
        Parameters      % struct — merged from returnDefaultParams + user overrides
        Activity        % (1 x N_cultures) struct array: .I .E .total .ratio .x .I_raw .E_raw .binned_mat
                        %   NOTE: .I/.E/.total are in spikes/bin (NOT Hz); divide by BinSize to convert
        BurstState      % (1 x N_cultures) cell: logical burst-state vector (true = burst bin)
        BurstCutouts    % (1 x N_cultures) struct: .all .accepted .rejected .clusters .cluster_idx .silhouette .accepted_locs
        Correlations    % (1 x N_cultures) struct: .pop_corr .cell_id
        NormalizedCutouts  % struct: .bursts (N_cultures x N_bins), .inh_frac (N_cultures x N_bins), .burst_counts
    end

    methods
        function eia = EIAnalyzer(ctc, parameters)
            arguments
                ctc        CellTypeClassifier
                parameters struct = struct()
            end
            assert(~isempty(ctc.UnitLabels), ...
                'CellTypeClassifier must be fully classified (run classifyUnits) before constructing EIAnalyzer');
            eia.Classifier  = ctc;
            eia.Parameters  = parseStructParameters(eia.returnDefaultParams(), parameters);
        end
    end

    methods (Access = private)
        function map = getCultureUnitMapping(eia)
        % GETCULTUREUNITMAPPING  Build culture-to-unit index map from FeatureStore metadata.
        %
        % Returns a struct with:
        %   .culture_ids     — string vector (N_units x 1) mapping each FS row to a culture ID
        %   .unique_cultures — ordered unique culture IDs
        %   .ud_order        — index into UnitDataArray for each FS row (0 = not found)
        %   .labels          — UnitLabels vector in FS row order
            ctc = eia.Classifier;
            fs  = ctc.FeatureStore;
            ud  = ctc.UnitDataArray;

            unit_rec_ids    = fs.UnitTable.RecordingID;
            culture_ids     = FeatureStore.getCultureIDsForUnits(unit_rec_ids, fs.MetadataTable, ctc.Parameters.CultureKeys);
            unique_cultures = unique(culture_ids, 'stable');

            unit_ids_fs   = fs.UnitTable.UnitID;
            ud_ids        = string({ud.UnitID})';
            [~, ud_order] = ismember(unit_ids_fs, ud_ids);

            map = struct( ...
                'culture_ids',     culture_ids, ...
                'unique_cultures', unique_cultures, ...
                'ud_order',        ud_order, ...
                'labels',          ctc.UnitLabels);
        end
    end

    methods (Static)
        function defaultParams = returnDefaultParams()
            defaultParams.Activity.BinSize          = 0.01;   % seconds per bin
            defaultParams.Activity.SecCutout        = [0, 1200];  % analysis window (s)
            defaultParams.Activity.SmoothWindow     = 5;      % movmean window (bins)
            defaultParams.Activity.RateType         = "mean"; % "mean" (per-unit average) | "sum" (total spikes)

            defaultParams.BurstDetection.MaxIERatio         = 2;    % burst when I/E ratio < this
            defaultParams.BurstDetection.NoiseFloorPercentile = 5;  % percentile of positive E values as noise floor
            defaultParams.BurstDetection.SmoothWindow       = 9;    % median filter width (bins)
            defaultParams.BurstDetection.PeakCutout         = 150;  % bins +/- burst peak to extract
            defaultParams.BurstDetection.PeakHeight         = 0.1;  % minimum peak height
            defaultParams.BurstDetection.PeakProminence     = 0.1;  % minimum peak prominence
            defaultParams.BurstDetection.MinPeakDistance    = 6;    % minimum bins between detected burst peaks
            defaultParams.BurstDetection.MaskPeaks          = false; % true: peak detection only in burst-state bins
            defaultParams.BurstDetection.MinBursts          = 4;    % minimum bursts required to process a culture
            defaultParams.BurstDetection.MaxAlignShift      = 50;   % maximum xcorr shift (bins) during alignment

            % Outlier rejection (Stage 1)
            defaultParams.BurstDetection.CorrelationThreshold = 0.8;  % min Pearson r with mean template
            defaultParams.BurstDetection.MaxFilterIter        = 5;    % max refinement iterations

            % Hierarchical clustering (Stage 2)
            defaultParams.BurstDetection.MaxClusters        = 4;    % maximum number of burst clusters to consider
            defaultParams.BurstDetection.MinSilhouette      = 0.3;  % min mean silhouette to accept k>1
            defaultParams.BurstDetection.ClusterDistance    = "correlation";  % pdist metric for clustering

            % Normalisation
            defaultParams.BurstDetection.NormalizationWindow = "full";  % "full" | "peak_to_end" (post-peak only)
            defaultParams.BurstDetection.TemporalFeatures    = false;  % append temporal context to burst clustering

            % Plotting
            defaultParams.Plotting.FontSize = 7;
        end
    end
end
