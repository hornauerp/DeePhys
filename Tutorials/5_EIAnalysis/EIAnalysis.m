%% Tutorial 5 — E/I Analysis
%
% Computes excitatory/inhibitory network activity traces, detects burst
% episodes, extracts aligned burst cutouts, and quantifies unit-population
% correlations from a trained CellTypeClassifier.
%
% Prerequisites:
%   - Tutorial 4 completed (CellTypeClassifier trained and saved)
%   - DeePhys on the MATLAB path

%% 1  Load a trained CellTypeClassifier
%
% If you saved the ctc via:
%   save('/path/to/ctc.mat', 'ctc');
% reload it here.  EIAnalyzer only requires ctc.UnitLabels to be set.

load('/path/to/ctc.mat', 'ctc');

% Verify the classifier is fully trained
assert(~isempty(ctc.UnitLabels), 'Run classifyUnits() before constructing EIAnalyzer.');

%% 2  Construct EIAnalyzer

eia = EIAnalyzer(ctc);

% With custom parameters:
%   params.Activity.BinSize          = 0.01;   % 10 ms bins
%   params.Activity.SecCutout        = [0, 600]; % first 600 s
%   params.BurstDetection.MaxIERatio = 2;       % I/E ratio to enter burst state
%   params.BurstDetection.SmoothWindow = 9;     % median filter width (bins)
%   eia = EIAnalyzer(ctc, params);

disp(eia.Parameters.Activity);
disp(eia.Parameters.BurstDetection);

%% 3  Compute per-culture activity traces
%
% For each culture, computes:
%   .E     — excitatory mean spike count per bin (spikes/bin; divide by BinSize for Hz)
%   .I     — inhibitory mean spike count per bin
%   .total — total classified-unit mean spike count per bin
%   .ratio — smoothed I/E ratio trace
%   .x     — time axis (s)

eia.computeActivity();

% Inspect the first culture
act1 = eia.Activity(1);
fprintf('Culture 1: %d time bins, E max %.4f spikes/bin, I max %.4f spikes/bin\n', ...
    numel(act1.x), max(act1.E), max(act1.I));

%% 4  Detect bursts
%
% Labels each time bin as burst (true) or non-burst (false). A bin is a burst
% when the I/E ratio is below MaxIERatio AND the excitatory rate exceeds the
% percentile-based noise floor.

eia.detectBursts();

% Inspect burst state for first culture
bs1 = eia.BurstState{1};
burst_frac = mean(bs1);  % BurstState is logical: true = burst
fprintf('Culture 1 burst fraction: %.1f%%\n', 100 * burst_frac);

%% 5  Extract aligned burst cutouts
%
% Detects burst peak times, extracts ±PeakCutout bins around each,
% and clusters them (t-SNE or correlation-based) to select canonical shapes.

eia.extractBurstCutouts();

bc1 = eia.BurstCutouts(1);
fprintf('Culture 1 bursts detected: %d\n', size(bc1.all.total, 1));

%% 6  Compute unit-population correlations
%
% For each unit, correlates its binned spike train with the population
% activity trace, separately for E and I populations.

eia.computeCorrelations();

corr1 = eia.Correlations(1);
fprintf('Culture 1 units with correlations: %d\n', numel(corr1.cell_id));

%% 7  Normalize burst cutouts across cultures

eia.normalizeBurstCutouts();

%% 8  Plotting — network activity

% Plot firing rate trace + burst state for culture 1
eia.PlotNetworkActivity(1);

%% 9  Plotting — E/I traces

% Plot stacked E and I firing rate traces for culture 1 (full recording)
eia.PlotEITraces(1);

% To zoom into a specific window, pass bin indices via StackCutout:
%   eia.PlotEITraces(1, 'StackCutout', 1:round(600/eia.Parameters.Activity.BinSize));

%% 10  Plotting — burst cutouts

eia.PlotBurstCutouts(1);

%% 10b  Plotting — burst-aligned raster
%
% Shows individual unit spike times aligned to burst peaks, sorted by cell type
% (excitatory first, then inhibitory) and coloured accordingly.

eia.PlotBurstRaster(1);

% To show only a specific burst or cluster:
%   eia.PlotBurstRaster(1, 'BurstIndex', 3);
%   eia.PlotBurstRaster(1, 'ClusterIndex', 1);

%% 11  Parameter tuning — burst detection
%
% Lower MaxIERatio = more bins flagged as bursts (excitatory-dominant threshold).
% Higher SmoothWindow = smoother ratio trace, less sensitive to short transients.

eia2 = EIAnalyzer(ctc);
eia2.Parameters.BurstDetection.MaxIERatio   = 1.5;
eia2.Parameters.BurstDetection.SmoothWindow = 15;
eia2.computeActivity();
eia2.detectBursts();
eia2.extractBurstCutouts();

%% 12  Decoupled workflow — EI analysis from saved classifier only
%
% You do not need to re-run the full classification pipeline to run EI
% analysis.  Load the saved ctc, verify UnitLabels, and proceed.

load('/path/to/ctc.mat', 'ctc');
eia3 = EIAnalyzer(ctc);
eia3.computeActivity();
eia3.detectBursts();
eia3.computeCorrelations();

%% 13  Multi-culture loop

n_cultures = numel(eia.Activity);
burst_fracs = zeros(1, n_cultures);
for c = 1:n_cultures
    if ~isempty(eia.BurstState{c})
        burst_fracs(c) = mean(eia.BurstState{c});  % logical: true = burst
    end
end
figure;
bar(burst_fracs);
xlabel('Culture index'); ylabel('Burst fraction');
title('Burst fraction per culture');
