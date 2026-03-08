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
    end

    methods (Static)
        function defaultParams = returnDefaultParams()
            defaultParams.Bootstrap.PreCutout   = [0, 1200];    % seconds — baseline window
            defaultParams.Bootstrap.PostCutout  = [6000, 7200]; % seconds — post-treatment window
            defaultParams.Bootstrap.BinSize     = 20;           % seconds per bin
            defaultParams.Bootstrap.NIter       = 1000;         % bootstrap iterations
            defaultParams.Bootstrap.Alpha       = 1e-10;        % significance level

            defaultParams.UMAP.NDims            = 2;            % UMAP output dimensions
            defaultParams.UMAP.NNeighbors       = 100;
            defaultParams.UMAP.MinDist          = 0.1;
            defaultParams.UMAP.Spread           = 1;
            defaultParams.UMAP.UnitFeatures     = ["FullACG","ReferenceWaveform"]; % feature groups to use
            defaultParams.UMAP.NormalizationVar = "ChipID";     % metadata field for within-group normalisation
            defaultParams.UMAP.GroupingVar      = "Concentration"; % metadata field grouping recordings (e.g. dose)
            defaultParams.UMAP.GroupingValues   = 0;            % value(s) of GroupingVar to use (0 = baseline)
            defaultParams.UMAP.TemplateDir      = tempdir;      % directory for UMAP template .mat file

            defaultParams.OutlierDetection.ContaminationFraction = 0.5;
            defaultParams.OutlierDetection.NObsPerLearner        = 50;
            defaultParams.OutlierDetection.DistancePercentile    = 95;
        end
    end
end
