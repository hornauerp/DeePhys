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
% Burst-event analysis (bursts from detectNetworkBursts, so the labels play no part in finding
% them; every metric per recording against random "I" sets matched in firing rate):
%   eia.computeNetworkBursts();   % onset / peak / offset / superbursts per recording
%   eia.computeBurstEI();         % E / I epoch rates, aligned traces, metrics vs null
%   eia.computeUnitTiming();      % per unit: latency to onset, participation, pre-onset rate
%   eia.computeBurstRhythm();     % burst-locked sub-peaks and shared fast fluctuations
%   T = eia.summarizeBurstEI();   % one row per recording: metadata, metrics, null, unit timing
%   eia.PlotBurstEI(r)            % aligned E / I traces and I/E with the null band
% The legacy detectBursts labels burst bins by I/E < MaxIERatio, so statements about the I/E
% ratio during those bursts are circular.
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
        NetworkBursts   % (1 x N_recordings) struct: RecordingID, Culture, UnitRows (FS rows), Duration, Bursts, Info
        BurstEI         % (1 x N_recordings) struct: RecordingID, Culture, Labels (of UnitRows), Result (burstEIMetrics)
        UnitTiming      % table: RecordingID, Culture, UnitRow, UnitID, Label + burstUnitTiming columns
        BurstRhythm     % (1 x N_recordings) struct: RecordingID, Culture, Result (burstRhythm)
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
        % GETCULTUREUNITMAPPING  Group FeatureStore rows into simultaneously recorded populations.
        %
        % Population activity is only meaningful for units recorded together, so groups are
        % recordings (RecordingID). A culture recorded several times (e.g. a dose series)
        % yields one group per recording; with baseline-only data a group is a culture.
        %
        % Returns a struct with:
        %   .culture_ids     — string vector (N_units x 1): group (recording) of each FS row
        %   .unique_cultures — ordered unique groups (RecordingIDs)
        %   .culture_names   — culture ID (CultureKeys) of each group, for labelling
        %   .ud_order        — index into UnitDataArray for each FS row (0 = not found)
        %   .labels          — labels in FS row order (Parameters.Labels; default UnitLabels)
            ctc = eia.Classifier;
            fs  = ctc.FeatureStore;
            ud  = ctc.UnitDataArray;

            unit_rec_ids    = string(fs.UnitTable.RecordingID);
            unique_recs     = unique(unit_rec_ids, 'stable');
            culture_of_unit = FeatureStore.getCultureIDsForUnits(unit_rec_ids, fs.MetadataTable, ctc.Parameters.CultureKeys);
            [~, first_row]  = ismember(unique_recs, unit_rec_ids);

            % UnitIDs repeat across recordings of one culture (and collide occasionally
            % across sortings): match UnitData by (UnitID, RecordingID).
            fs_keys       = string(fs.UnitTable.UnitID) + "|" + unit_rec_ids;
            ud_keys       = string({ud.UnitID})' + "|" + string({ud.RecordingID})';
            [~, ud_order] = ismember(fs_keys, ud_keys);

            map = struct( ...
                'culture_ids',     unit_rec_ids, ...
                'unique_cultures', unique_recs, ...
                'culture_names',   culture_of_unit(first_row), ...
                'ud_order',        ud_order, ...
                'labels',          eia.unitLabels());
        end

        function labels = unitLabels(eia)
        % UNITLABELS  Labels of the E / I populations (Parameters.Labels), FeatureStore row order.
            pl = eia.Parameters.Labels;
            ctc = eia.Classifier;
            switch lower(string(pl.Source))
                case "graph"
                    labels = ctc.UnitLabels;
                case "composition"
                    assert(isfield(ctc.Composition, 'Labels'), 'EIAnalyzer:noComposition', 'Run ctc.estimateComposition() first.');
                    labels = ctc.Composition.Labels;
                case "confident"
                    [~, labels] = evalc('ctc.confidentLabels(MinProbability=pl.MinProbability, MinSpikes=pl.MinSpikes)');
                case "responders"
                    labels = eia.responderLabels();
                otherwise
                    error('EIAnalyzer:labelSource', 'Labels.Source must be "graph", "composition", "confident" or "responders".');
            end
            labels = labels(:)';
        end

        function labels = responderLabels(eia)
        % RESPONDERLABELS  Ground-truth populations from the perturbation (the paper's Fig 4A-C):
        % I = units whose ResponsiveUnitDirection is Labels.ResponderDirection, E = every other unit
        % of the same cultures. Cultures: those matching Labels.ResponderFilter ({field, value},
        % e.g. {'Construct', 1}) or, if empty, those with at least one responder; NaN elsewhere.
            pl = eia.Parameters.Labels; ctc = eia.Classifier; fs = ctc.FeatureStore;
            assert(~isempty(ctc.ResponsiveUnitDirection), 'EIAnalyzer:noResponders', ...
                'Run ctc.identifyResponsiveUnits() first (Labels.Source "responders").');
            dir_ = string(ctc.ResponsiveUnitDirection(:)');
            rec = string(fs.UnitTable.RecordingID)';
            cult = FeatureStore.getCultureIDsForUnits(rec(:), fs.MetadataTable, ctc.Parameters.CultureKeys)';
            is_i = dir_ == string(pl.ResponderDirection);
            if isempty(pl.ResponderFilter)
                in_set = ismember(cult, unique(cult(is_i)));
            else
                meta = fs.MetadataTable; v = meta.(pl.ResponderFilter{1});
                if iscell(v), v = string(v); end
                in_set = ismember(rec, string(meta.RecordingID(v == pl.ResponderFilter{2})));
                in_set = ismember(cult, unique(cult(in_set)));      % every recording of those cultures
            end
            labels = nan(1, numel(dir_)); labels(in_set) = 1; labels(in_set & is_i) = 2;
        end
    end

    methods (Static)
        function defaultParams = returnDefaultParams()
            % Labels of the E and I populations:
            %   "graph" (default)  — UnitLabels of classifyUnits (the published analysis)
            %   "composition"      — Composition.Labels (estimateComposition, culture-aware)
            %   "confident"        — confidentLabels: only units with P(inhibitory) >= MinProbability
            %                        or <= 1 - MinProbability (and >= MinSpikes spikes); the rest are
            %                        left out of both populations
            %   "responders"       — ground truth of a perturbation: I = units responding in
            %                        ResponderDirection, E = the other units of those cultures
            %                        (ResponderFilter {field, value}, or cultures with a responder)
            defaultParams.Labels.Source         = "graph";
            defaultParams.Labels.MinProbability = 0.8;
            defaultParams.Labels.MinSpikes      = 0;    % "confident": minimum spikes behind the features
            defaultParams.Labels.ResponderDirection = "increase";   % "responders": hM3Dq in interneurons
            defaultParams.Labels.ResponderFilter    = {};           % e.g. {'Construct', 1}

            defaultParams.Activity.BinSize          = 0.01;   % seconds per bin
            defaultParams.Activity.SecCutout        = [0, 1200];  % analysis window (s); ends at the recording's end if earlier (Inf = whole recording)
            defaultParams.Activity.SmoothWindow     = 5;      % movmean window (bins)
            defaultParams.Activity.RateType         = "mean"; % "mean" (per-unit average) | "sum" (total spikes)
            % RatioMethod: "smoothed" (default) = movmean(I) ./ movmean(E), NaN where the window
            %   holds no E spike; "perbin" = movmean of the per-bin I/E over bins with E > 0 (the
            %   published method; with sparse 10 ms counts it is biased by count noise: inflated
            %   where counts are high, deflated where they are low, even at a constant true I/E).
            defaultParams.Activity.RatioMethod      = "smoothed";

            % MaxIERatio: burst when I/E ratio < this (and E above the noise floor). In
            %   simulated bursting networks both RatioMethods separate burst bins equally well
            %   at the same threshold ("smoothed" slightly better); the E noise floor does most of
            %   the work, and thresholds near the in-burst I/E cut into burst bins.
            defaultParams.BurstDetection.MaxIERatio         = 2;
            defaultParams.BurstDetection.NoiseFloorPercentile = 5;  % percentile of positive E values as noise floor
            defaultParams.BurstDetection.SmoothWindow       = 9;    % median filter width (bins)
            defaultParams.BurstDetection.PeakCutout         = 150;  % bins +/- burst peak to extract
            defaultParams.BurstDetection.PeakHeight         = 0.1;  % minimum peak height (fraction of the activity maximum)
            defaultParams.BurstDetection.PeakProminence     = 0.1;  % minimum peak prominence (fraction of the activity maximum)
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

            % Burst-event analysis
            defaultParams.NetworkBursts.Options = struct();   % name-value options of detectNetworkBursts
            defaultParams.BurstEI.PreWindow      = 0.3;       % s before onset
            defaultParams.BurstEI.PostWindow     = 0.3;       % s after offset
            defaultParams.BurstEI.Margin         = 0.5;       % interburst reference: outside bursts +- Margin
            defaultParams.BurstEI.MinBetween     = 10;        % s of interburst time needed for the folds
            defaultParams.BurstEI.BinSize        = 0.01;      % aligned traces (s)
            defaultParams.BurstEI.Window         = [-1.5 1.5];
            defaultParams.BurstEI.NNull          = 200;       % rate-matched random "I" sets
            defaultParams.BurstEI.NullStrata     = 5;         % firing-rate quantile strata
            defaultParams.BurstEI.MinUnits       = 3;         % per population
            defaultParams.BurstEI.MinBursts      = 5;
            defaultParams.BurstEI.Seed           = 1;
            defaultParams.UnitTiming.PreWindow   = 0.3;
            defaultParams.UnitTiming.PostWindow  = 0.3;
            defaultParams.BurstRhythm.BinSize      = 0.002;
            defaultParams.BurstRhythm.SmoothSigma  = 0.025;
            defaultParams.BurstRhythm.LockedWindow = [-0.1 0.3];
            defaultParams.BurstRhythm.SharedWindow = 0.25;
            defaultParams.BurstRhythm.NSplits      = 5;
            defaultParams.BurstRhythm.NShift       = 20;
            defaultParams.BurstRhythm.MinBursts    = 10;

            % Plotting
            defaultParams.Plotting.FontSize = 7;
        end
    end
end
