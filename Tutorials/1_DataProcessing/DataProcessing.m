%% Tutorial 1 — Data Processing Pipeline
%
% Covers the full pipeline from raw Kilosort output to a FeatureStore
% ready for analysis: SpikeData → RecordingProcessor → FeatureStore.
%
% Prerequisites: DeePhys on the MATLAB path (run startup.m from repo root).
%
% Fill in the placeholder paths in Section 1 before running.

%% 1  Setup — fill in your paths

% Path to one Kilosort output directory (must contain spike_times.npy etc.)
ks_path = '/path/to/kilosort_output';

% Optional: path to a parent (concatenated) recording.
% Set to '' if your recording was not split from a longer session.
parent_ks_path = '/path/to/parent_kilosort_output';

% Where to save the processed RecordingProcessor and FeatureStore files
save_dir = '/path/to/save';

%% 2  Metadata struct
%
% One struct per recording. Any scalar fields are stored in the FeatureStore
% and become available for subsetting and ML labels.

metadata = struct();
metadata.ChipID       = 'Chip001';
metadata.PlatingDate  = '2024-01-15';
metadata.RecordingDate = '2024-02-05';  % DIV computed automatically if absent
metadata.DIV          = 21;
metadata.Mutation     = 'WT';
metadata.Concentration = 0;

% If split from a parent recording, specify the parent path here.
% If omitted, SpikeData tries to auto-detect a sibling 'qc_output' folder.
if ~isempty(parent_ks_path)
    metadata.ParentInputPath = parent_ks_path;
end

%% 3  Load Kilosort output into SpikeData

sd = SpikeData.fromKilosort(ks_path, metadata);

% Inspect
fprintf('RecordingID : %s\n', sd.RecordingID);
fprintf('Duration    : %.1f s\n', sd.Duration);
fprintf('N spikes    : %d\n', numel(sd.SpikeTimes));
fprintf('N templates : %d\n', size(sd.TemplateWaveforms, 1));
fprintf('Parent path : %s\n', sd.ParentPath);

%% 4  Create RecordingProcessor (no analysis yet)
%
% Pass optional parameter overrides in the second argument (struct).
% All unspecified fields fall back to RecordingProcessor.returnDefaultParams().

proc = RecordingProcessor(sd);

% Alternatively, load directly from disk in one step:
%   proc = RecordingProcessor.fromKilosort(ks_path, metadata);

%% 5  Quality control — filter templates into proc.Units

proc.runQC();

fprintf('Units after QC: %d\n', numel(proc.Units));

% Inspect individual units
if ~isempty(proc.Units)
    u = proc.Units(1);
    fprintf('  Unit 1 — TemplateID: %d, N spikes: %d, FR: %.2f Hz\n', ...
        u.TemplateID, numel(u.SpikeTimes), numel(u.SpikeTimes) / sd.Duration);
end

%% 6  Compute per-unit features
%
% Computes ActivityFeatures, WaveformFeatures, RegularityFeatures, Catch22.
% Results stored in proc.UnitFeatureTable.

proc.computeUnitFeatures();

% Peek at the table
disp(proc.UnitFeatureTable(1:min(5, height(proc.UnitFeatureTable)), :));

%% 7  Compute parent features (full-recording ACGs)
%
% Reads spike data from proc.SpikeData.ParentPath, computes CCG once for ALL
% templates in the parent recording (shared cache across siblings), then
% stores the autocorrelogram diagonal as Parent_ACG1..N columns.
%
% If no parent path is available, a warning is printed and the step is skipped
% silently — child ACG features serve as automatic fallback.

proc.computeParentFeatures();

% Parent ACG columns are now in proc.UnitFeatureTable with 'Parent_' prefix
parent_cols = startsWith(string(proc.UnitFeatureTable.Properties.VariableNames), 'Parent_');
fprintf('Parent feature columns added: %d\n', sum(parent_cols));

%% 8  Compute network features

proc.computeNetworkFeatures();

disp(proc.NetworkFeatureTable);

%% 9  Compute connectivity (CCG / STTC)
%
% Monosynaptic connection windows are now configurable via the CCG sub-struct
% of proc.Parameters.Connectivity. The defaults match the values used in
% previous DeePhys releases; override only when needed for a different
% preparation (e.g. iPSC cultures vs. acute slices):
%
%   proc.Parameters.Connectivity.CCG.PostStart  = 0.0008;  % post-spike excit. window start (s) — default 0.8 ms
%   proc.Parameters.Connectivity.CCG.PostEnd    = 0.004;   % post-spike excit. window end   (s) — default 4.0 ms
%   proc.Parameters.Connectivity.CCG.PreStart   = 0.0032;  % pre-spike baseline window start (s) — default 3.2 ms
%   proc.Parameters.Connectivity.CCG.BonfWindow = 0.005;   % Bonferroni-correction window    (s) — default 5.0 ms

proc.computeConnectivity();

%% 10  Compute cell-type features (E/I-stratified)
%
% Produces E/I-stratified graph, balance, activity, burst, and correlation
% features — columns with suffixes like _EE, _II, _EI, or prefixes like
% ExcitatoryFraction, MeanFiringRate_E/I.
%
% Requires CellTypeLabels to be set (populated by Tutorial 4 after running
% CellTypeClassifier.classifyUnits and RecordingProcessor.applyLabelsFromClassifier).
% If CellTypeLabels is empty or all-NaN, this step returns immediately with no
% new columns added.
%
% Feature groups produced (accessible via fs.recordingMatrix("CellTypeBalance") etc.):
%   CellTypeGraph       — graph metrics on E-only and I-only sub-networks
%   CellTypeBalance     — ExcitatoryFraction, MeanFiringRate_E/I, FiringRateRatio_EI
%   CellTypeActivity    — MeanCV2_E/I, MeanFanoFactor_E/I, MeanLvR_E/I
%   CellTypeBurst       — burst lead/participation per cell type
%   CellTypeCorrelation — MeanSTTC_EE/II/EI, STTCSynchronyIndex

proc.computeCellTypeFeatures();
fprintf('CellTypeFeatures status: %s\n', proc.Status.CellTypeFeatures);

%% 11  Compute spatial analysis
%
% Computes spatial distribution features from electrode coordinates stored in
% proc.SpikeData.ElectrodeCoordinates. Skipped with a warning if coordinates
% are absent. Enriched by CellTypeLabels when available.
%
% Sub-analyses run automatically:
%   A. Spatial spread    — ConvexHullArea, ChipCoverage, MeanPairwiseDistance, CentroidSpread
%   B. E/I mixing        — SpatialMixingIndex, MeanFractionExcNN_E/I (requires labels)
%   D. FC decay          — SpatialDecayTau, DistanceFCCorrelation (requires Connectivity)
%   E. FR gradient       — SpatialFRMoransI, CenterPeripheryFR_ratio
%   F. E/I balance map   — SpatialEIVariability, SpatialEIMoransI (requires labels)
%   G. Ripley's K/L      — RipleysL_max, ClusterScale (optional)
%   H. Burst propagation — BurstOriginDispersion, MeanBurstPropagationSpeed (requires Bursts)
%
% All spatial columns are accessible via fs.recordingMatrix("SpatialFeatures")
% or fs.unitMatrix("SpatialFeatures") after FeatureStore assembly.

proc.computeSpatialAnalysis();
fprintf('SpatialFeatures status: %s\n', proc.Status.SpatialFeatures);

% Inspect spatial columns added to the network feature table
if ~isempty(proc.NetworkFeatureTable)
    sp_names = string(proc.NetworkFeatureTable.Properties.VariableNames);
    sp_prefixes = ["ConvexHull","ChipCoverage","MeanPairwise","CentroidSpread", ...
                   "SpatialMixing","MeanFractionExc","SpatialFR","CenterPeriphery", ...
                   "SpatialDecay","DistanceFC","SpatialEI","BurstOrigin", ...
                   "MeanBurstProp","Ripley","ClusterScale"];
    sp_mask = false(size(sp_names));
    for p = sp_prefixes, sp_mask = sp_mask | startsWith(sp_names, p); end
    fprintf('Spatial network features added: %d columns\n', sum(sp_mask));
end

% --- Visualization ---
%
% Two standalone functions are available for spatial inspection:
%
%   plotSpatialUnitMap(proc.SpikeData.ElectrodeCoordinates, proc.Units, ...
%       proc.CellTypeLabels);
%   % → Scatter map of units at reference electrode positions, colored by cell type.
%   %   Pass [] as third argument to skip cell-type coloring.
%
%   [~, viz_data] = computeSpatialEIBalance( ...
%       proc.SpikeData.ElectrodeCoordinates, proc.Units, proc.CellTypeLabels);
%   plotSpatialEIBalance(viz_data, proc.SpikeData.ElectrodeCoordinates);
%   % → Kernel-smoothed heatmap of local excitatory fraction across the MEA chip.

%% 12  Shortcut: run all steps in sequence

proc2 = RecordingProcessor(sd);
proc2.runAll();   % runQC + computeUnitFeatures + computeParentFeatures
                  %       + computeNetworkFeatures + computeConnectivity
                  %       + computeCellTypeFeatures + computeSpatialAnalysis

%% 13  Save and load

proc_file = fullfile(save_dir, 'RecordingProcessor.mat');
proc.save(proc_file);

proc_loaded = RecordingProcessor.load(proc_file);
fprintf('Loaded %d units from disk.\n', numel(proc_loaded.Units));

% Inspect status flags
disp(proc_loaded.Status);

%% 14  Batch loading with parallel workers

% Supply a cell array or string array of saved RecordingProcessor .mat paths
proc_paths = {proc_file};   % replace with your full list
procs = RecordingProcessor.loadMany(proc_paths);
fprintf('Batch-loaded %d processors.\n', numel(procs));

%% 15  Assemble FeatureStore from multiple processors

fs = FeatureStore.fromProcessors(procs);

% Three tables
fprintf('UnitTable     : %d rows × %d cols\n', height(fs.UnitTable),      width(fs.UnitTable));
fprintf('RecordingTable: %d rows × %d cols\n', height(fs.RecordingTable), width(fs.RecordingTable));
fprintf('MetadataTable : %d rows × %d cols\n', height(fs.MetadataTable),  width(fs.MetadataTable));

%% 16  Inspect table structure

% Column names in UnitTable
disp(string(fs.UnitTable.Properties.VariableNames)');

% First few rows
disp(fs.UnitTable(1:min(5, height(fs.UnitTable)), 1:8));

%% 17  Save and load FeatureStore

fs_file = fullfile(save_dir, 'FeatureStore.mat');
fs.save(fs_file);

fs2 = FeatureStore.load(fs_file);
fprintf('Loaded FeatureStore: %d units, %d recordings.\n', ...
    height(fs2.UnitTable), height(fs2.RecordingTable));
