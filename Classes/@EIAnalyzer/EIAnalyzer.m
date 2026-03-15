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
%   eia.extractBurstCutouts();    % aligned burst cutouts, t-SNE clustered
%   eia.computeCorrelations();    % unit ↔ population E/I correlations
%
% Plotting:
%   eia.PlotNetworkActivity(c)    % firing rate coloured by burst state
%   eia.PlotBurstCutouts(c)       % mean E/I burst shapes and gradients
%   eia.PlotEITraces(c, range)    % stacked E/I + total firing rate

    properties
        Classifier      % CellTypeClassifier — provides UnitLabels, UnitList, RecordingGroup
        Parameters      % struct — merged from returnDefaultParams + user overrides
        Activity        % (1 x N_cultures) struct array: .I .E .total .ratio .x
        BurstState      % (1 x N_cultures) cell: binary burst-state vector per culture
        BurstCutouts    % (1 x N_cultures) struct: .total .i_mat .e_mat (selected variant)
        Correlations       % (1 x N_cultures) struct: .pop_corr .cell_id
        NormalizedCutouts  % struct: .bursts (N_cultures x N_bins), .ie (N_cultures x N_bins)
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

    methods (Static)
        function defaultParams = returnDefaultParams()
            defaultParams.Activity.BinSize          = 0.01;       % seconds per bin
            defaultParams.Activity.SecCutout        = [0, 1200];  % analysis window (s)

            defaultParams.BurstDetection.Threshold       = 2;     % I/E ratio threshold for burst state
            defaultParams.BurstDetection.SmoothWindow    = 9;     % median filter width (bins)
            defaultParams.BurstDetection.PeakCutout      = 150;   % bins ± burst peak to extract
            defaultParams.BurstDetection.PeakHeight      = 0.1;   % minimum peak height
            defaultParams.BurstDetection.PeakProminence  = 0.1;   % minimum peak prominence
            defaultParams.BurstDetection.MinPeakDistance = 6;     % minimum bins between detected burst peaks
            defaultParams.BurstDetection.MaskPeaks       = false; % true: peak detection only in burst-state bins (z==2)
                                                                   % false: peak detection on full activity trace (main-branch default)
            defaultParams.BurstDetection.SelectedVariant = 4;     % burst variant: 1=all, 2=cluster1, 3=cluster2, 4=both-aligned
            defaultParams.BurstDetection.BurstFiltering  = "correlation"; % "tsne" | "correlation" | "none"
            defaultParams.BurstDetection.CorrelationThreshold = 0.8;  % min Pearson r with mean template (for "correlation" filtering)
            defaultParams.BurstDetection.MaxFilterIter        = 5;    % max refinement iterations (for "correlation" filtering)
            defaultParams.BurstDetection.RNGSeed              = 42;   % seed for t-SNE + kmedoids reproducibility
        end
    end
end
