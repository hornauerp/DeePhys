classdef CellTypeClassifier < handle
% CELLTYPECLASSIFIER  Supervised cell-type classification pipeline for MEA perturbation experiments.
%
% Identifies inhibitory (interneuron) vs excitatory neurons in MEA cultures
% using a bootstrap firing-rate response test followed by supervised UMAP classification.
%
% USAGE:
%   ctc = CellTypeClassifier(rg, params);
%   ctc.identifyResponsiveUnits();    % bootstrap pre/post firing rate comparison
%   ctc.generateTrainLabels();        % UMAP embedding + isolation forest label generation
%   ctc.classifyUnits();              % supervised UMAP projection and classification
%   labels = ctc.UnitLabels;          % 1 = excitatory, 2 = inhibitory, NaN = unclassified

    properties
        RecordingGroup          % RecordingGroup with .Cultures already set
        Parameters              % struct — merged from returnDefaultParams + user overrides
        UnitList                % (1 x N) Unit array — all units across all cultures
        ResponsiveUnitIdx       % (1 x N) logical — units identified as responsive (inhibitory candidates) by bootstrap
        TrainLabels             % struct: .sorted_train_ids, .sorted_y_train, .umap_train_idx, .umap_test_idx
        UnitLabels              % (1 x N) double — 1=excitatory, 2=inhibitory, NaN=unclassified
        UMAP                    % trained UMAP model returned by run_umap
        Reduction = struct('Unsupervised', [], 'Train', [], 'Test', [], 'External', [])
        HarmonizedWaveforms     % (N_samples × N_units) processed waveforms (upsampled, aligned, trimmed)
        HarmonizedACGs          % (N_bins × N_units)   ACGs at Harmonization params
        HarmonizedSR            % scalar — waveform sampling rate after harmonization
        CachedExtraction        % struct caching extractUnitWaveformsAndACGs output (see getOrExtract)
                                % struct storing UMAP embeddings:
                                %   .Unsupervised — (N x D) from generateTrainLabels (all units, unsupervised)
                                %   .Train        — (N_train x D) from classifyUnits (supervised training embedding)
                                %   .Test         — (N_test  x D) from classifyUnits (supervised test projection)
                                %   .External     — (N_ext  x D) from classifyExternalUnits
    end

    methods
        function ctc = CellTypeClassifier(rg, parameters)
            arguments
                rg RecordingGroup
                parameters struct = struct()
            end
            ctc.RecordingGroup = rg;
            ctc.Parameters = parseStructParameters(ctc.returnDefaultParams(), parameters);
            % UnitList is populated by identifyResponsiveUnits(), which also sets
            % the cumulative culture offsets needed for ResponsiveUnitIdx.
        end
    

        function [wf, acg, sr] = getOrExtract(ctc, unit_list)
        % getOrExtract  Return cached extraction or compute and cache.
        %   Avoids redundant FullACG recomputation across generateTrainLabels,
        %   classifyUnits, etc. Cache is invalidated when unit count or
        %   harmonization parameters change.
            ph = ctc.Parameters.Harmonization;
            cache = ctc.CachedExtraction;
            if ~isempty(cache) ...
                    && cache.N == numel(unit_list) ...
                    && abs(cache.ACGBinSize - ph.ACGBinSize) < 1e-12 ...
                    && abs(cache.ACGLag - ph.ACGLag) < 1e-12 ...
                    && cache.WaveformTargetSR == ph.WaveformTargetSamplingRate ...
                    && cache.ACGSource == ph.ACGSource
                wf  = cache.wf;
                acg = cache.acg;
                sr  = cache.sr;
                fprintf('Using cached waveform/ACG extraction (%d units)\n', cache.N);
                return
            end
            [wf, acg, sr] = extractUnitWaveformsAndACGs(ctc, unit_list);
            ctc.CachedExtraction = struct( ...
                'wf', wf, 'acg', acg, 'sr', sr, ...
                'N', numel(unit_list), ...
                'ACGBinSize', ph.ACGBinSize, ...
                'ACGLag', ph.ACGLag, ...
                'WaveformTargetSR', ph.WaveformTargetSamplingRate, ...
                'ACGSource', ph.ACGSource);
        end
    end

    methods (Static)
        function defaultParams = returnDefaultParams()
            defaultParams.Bootstrap.PreCutout   = [0, 1200];    % seconds — baseline window
            defaultParams.Bootstrap.PostCutout  = [6000, 7200]; % seconds — post-treatment window
            defaultParams.Bootstrap.BinSize     = 20;           % seconds per bin
            defaultParams.Bootstrap.NIter       = 1000;         % bootstrap iterations
            defaultParams.Bootstrap.Alpha       = 1e-10;        % significance level

            defaultParams.UMAP.NDims            = 10;           % unsupervised UMAP output dimensions (label generation)
            defaultParams.UMAP.NNeighbors       = 50;           % unsupervised UMAP n_neighbors
            defaultParams.UMAP.MinDist          = 0.1;
            defaultParams.UMAP.Spread           = 1;
            defaultParams.UMAP.SupervisedNDims      = 2;        % supervised UMAP output dimensions (classification)
            defaultParams.UMAP.SupervisedNNeighbors = 100;      % supervised UMAP n_neighbors
            defaultParams.UMAP.UnitFeatures     = ["FullACG","ReferenceWaveform"]; % feature groups to use
            defaultParams.UMAP.NormalizationVar = "ChipID";     % metadata field for within-group normalisation
            defaultParams.UMAP.GroupingVar      = "Concentration"; % metadata field grouping recordings (e.g. dose)
            defaultParams.UMAP.GroupingValues   = 0;            % value(s) of GroupingVar to use (0 = baseline)
            defaultParams.UMAP.TrainingCultureIdx = [];         % culture indices for label generation UMAP (empty = all)
            defaultParams.UMAP.TemplateDir      = tempdir;      % directory for UMAP template .mat file

            defaultParams.OutlierDetection.ContaminationFraction = 0.5;
            defaultParams.OutlierDetection.NObsPerLearner        = 50;
            defaultParams.OutlierDetection.DistancePercentile    = 80;
            defaultParams.OutlierDetection.CounterexampleRatio   = 2;   % excitatory:inhibitory sampling ratio

            defaultParams.RNGSeed = 42;  % seed for reproducible label generation & classification

            % Harmonization parameters — shared target spec for DeePhys and external data.
            % buildFeatureMatrix and extractUnitWaveformsAndACGs both read these to
            % ensure DeePhys and external units are projected into the same feature space.
            % Change these to harmonize to a different ACG resolution or waveform rate.
            defaultParams.Harmonization.WaveformTargetSamplingRate = 300000;  % Hz — target rate after interpolation (must be integer multiple of recording SR)
            defaultParams.Harmonization.ACGBinSize                 = 0.0005; % s — ACG bin width (0.5 ms)
            defaultParams.Harmonization.ACGLag                     = 0.1;    % s — one-sided ACG lag (100 ms → 401 bins)
            defaultParams.Harmonization.ACGSource                  = "FullACG"; % "FullACG" or "ACG" — which DeePhys ACG property to use

            % Fallback parameters for fully external workflows (ctc.UnitList is empty).
            % When ctc.UnitList is non-empty, WaveformSamplingRate and WaveformNSamples
            % are inferred from the actual DeePhys units instead.
            defaultParams.External.WaveformSamplingRate = 30000;  % Hz — native recording rate
            defaultParams.External.WaveformNSamples     = 82;     % samples at WaveformSamplingRate (~2.73 ms)
        end
    end
end
