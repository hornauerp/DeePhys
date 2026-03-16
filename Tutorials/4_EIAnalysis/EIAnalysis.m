% EIAnalysis.m
%
% Excitatory/Inhibitory network analysis using @EIAnalyzer. Requires a
% trained CellTypeClassifier with unit labels. Computes per-culture E/I
% activity traces, burst detection via I/E ratio thresholding, aligned
% burst cutout extraction, and unit-population correlations.
%
% This class replaces the manual E/I analysis scripts that previously
% required assembling spike rasters, computing firing rate traces, and
% detecting bursts in separate steps.
%
% PREREQUISITES:
%   - A trained CellTypeClassifier (see 3_CellTypeClassification tutorial)
%   - run_umap toolbox on path

%% 0 - Setup
addpath(genpath("/path/to/DeePhys")) % CHANGE

%% 1 - Load a trained CellTypeClassifier
% Either run the CellTypeClassification tutorial first, or load a saved one:
% load('ctc_results.mat', 'ctc');

% For this tutorial, assume 'ctc' is already in the workspace with
% ctc.UnitLabels populated (1=excitatory, 2=inhibitory).
assert(~isempty(ctc.UnitLabels), 'CellTypeClassifier must be classified first')

%% 2 - Initialise EIAnalyzer
ei_params = struct();

% Activity trace parameters
ei_params.Activity.BinSize       = 0.01;       % 10 ms bins
ei_params.Activity.SecCutout     = [0, 1200];  % analysis window (s)

% Burst detection (based on I/E ratio threshold)
ei_params.BurstDetection.Threshold      = 2;    % I/E ratio above this = burst
ei_params.BurstDetection.SmoothWindow   = 9;    % median filter width
ei_params.BurstDetection.PeakCutout     = 150;  % bins around each burst peak
ei_params.BurstDetection.PeakHeight     = 0.1;
ei_params.BurstDetection.PeakProminence = 0.1;
ei_params.BurstDetection.BurstFiltering = "correlation"; % "tsne", "correlation", "none"
ei_params.BurstDetection.CorrelationThreshold = 0.8;

eia = EIAnalyzer(ctc, ei_params);

%% 3 - Compute E/I activity traces
% For each culture, computes excitatory, inhibitory, and total firing rate
% traces in the specified time window.
eia.computeActivity();

% Results stored in eia.Activity(c):
%   .E     — excitatory firing rate trace
%   .I     — inhibitory firing rate trace
%   .total — total firing rate trace
%   .ratio — I/E ratio trace
%   .x     — time vector

%% 4 - Visualise activity traces
% Plot for the first culture
c = 1;
figure('Color', 'w');
eia.PlotEITraces(c);

% Or with a specific time range
figure('Color', 'w');
eia.PlotEITraces(c, [0, 300]); % first 5 minutes

%% 5 - Detect bursts
% Binary burst-state labelling based on I/E ratio threshold.
% When I/E > threshold, the network is in a "burst state".
eia.detectBursts();

% Visualise with burst states colour-coded
figure('Color', 'w');
eia.PlotNetworkActivity(c);

%% 6 - Extract burst cutouts
% Aligns burst profiles to their peak and optionally filters by
% correlation with the mean template.
eia.extractBurstCutouts();

% Results in eia.BurstCutouts(c):
%   .total  — (N_bursts x N_bins) total activity around peak
%   .i_mat  — (N_bursts x N_bins) inhibitory activity
%   .e_mat  — (N_bursts x N_bins) excitatory activity

%% 7 - Visualise burst profiles
figure('Color', 'w');
eia.PlotBurstCutouts(c);
% Shows mean E/I burst shapes with SEM and their temporal derivatives

%% 8 - Unit-population correlations
% Computes correlation between each unit's activity and the population
% E or I trace to identify strongly coupled cells.
eia.computeCorrelations();

% Results in eia.Correlations(c):
%   .pop_corr — (N_units x 1) correlation with population trace
%   .cell_id  — unit indices

%% 9 - Cross-culture normalised burst profiles
% Align burst profiles across all cultures for group comparison.
eia.normalizeBurstCutouts();

% Results in eia.NormalizedCutouts:
%   .bursts — (N_cultures x N_bins) normalised burst traces
%   .ie     — (N_cultures x N_bins) normalised I/E traces

%% 10 - Summary
% The EIAnalyzer object stores all results and can be saved:
% save('eia_results.mat', 'eia');

% Key outputs for downstream analysis:
%   eia.Activity       — per-culture E/I/total firing rate traces
%   eia.BurstState     — binary burst-state vectors
%   eia.BurstCutouts   — aligned burst profiles (E, I, total)
%   eia.Correlations   — unit-population coupling
