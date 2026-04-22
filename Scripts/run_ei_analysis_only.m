% run_ei_analysis_only.m
%
% E/I network analysis from a pre-computed CellTypeClassifier.
%
% Use when classification is already done and only EIAnalyzer parameters
% need adjusting (burst threshold, analysis window, etc.).

%% 0 — Setup

deephys_root = getenv('DEEPHYS_ROOT');
if isempty(deephys_root)
    deephys_root = fileparts(fileparts(mfilename('fullpath')));
end
addpath(genpath(deephys_root));

%% 1 — Load pre-computed CellTypeClassifier

% Set save_path to the .mat file produced by run_cell_type_classification.m
save_path = '/path/to/cell_type_classification_results.mat';
loaded = load(save_path, 'ctc');
ctc    = loaded.ctc;

fprintf('Loaded: %d units (%d exc, %d inh)\n', ...
    numel(ctc.UnitLabels), sum(ctc.UnitLabels == 1), sum(ctc.UnitLabels == 2));

%% 2 — Configure EIAnalyzer

eia_params = struct();
eia_params.Activity.BinSize              = 0.01;
eia_params.Activity.SecCutout            = [0, 1200];
eia_params.BurstDetection.Threshold      = 2;
eia_params.BurstDetection.SmoothWindow   = 9;
eia_params.BurstDetection.PeakCutout     = 150;
eia_params.BurstDetection.PeakHeight     = 0.1;
eia_params.BurstDetection.PeakProminence = 0.1;
eia_params.BurstDetection.MinPeakDistance = 6;
eia_params.BurstDetection.SelectedVariant = 4;

eia = EIAnalyzer(ctc, eia_params);

%% 3 — Run analysis

eia.computeActivity();
fprintf('Activity computed for %d cultures.\n', numel(eia.Activity));

eia.detectBursts();
eia.extractBurstCutouts();
eia.computeCorrelations();
eia.normalizeBurstCutouts();

%% 4 — Visualise

n_cultures  = numel(eia.Activity);
norm_bursts = eia.NormalizedCutouts.bursts;
norm_ie     = eia.NormalizedCutouts.ie;

for c = 1:n_cultures
    eia.PlotNetworkActivity(c);
    drawnow;
end

figure('Color', 'w');
tiledlayout(2, 3, 'TileSpacing', 'compact');

nexttile; imagesc(norm_ie);     title('I/E ratio — all cultures',      'FontWeight', 'normal'); colorbar; axis tight
nexttile; imagesc(norm_bursts); title('Total activity — all cultures', 'FontWeight', 'normal'); colorbar; axis tight
nexttile;
plot(mean(norm_ie, 'omitnan'), 'r'); hold on;
plot(mean(norm_bursts, 'omitnan'), 'k');
legend('I/E', 'Total', 'Location', 'northwest'); xlabel('Bin (10 ms)'); box off

%% 5 — Save

save_dir_eia = fileparts(save_path);
if ~isfolder(save_dir_eia); mkdir(save_dir_eia); end
save_path_eia = fullfile(save_dir_eia, 'ei_analysis_results');
save(save_path_eia, 'eia');
fprintf('EIAnalyzer saved to: %s.mat\n', save_path_eia);
