% FeatureInspection.m
%
% After feature extraction, inspect individual MEArecording objects to
% verify QC, examine connectivity, waveforms, and burst detection results.
%
% WHAT CHANGED (2026 refactor):
%   - Unit features now include CV2, LV, LvR, FanoFactor, HalfWidth.
%   - Graph features include Modularity and SmallWorldness.
%   - Burst detection supports 'logISI' method alongside the default 'ISIN'.
%   - Zero-imputation replaced with NaN + median imputation. Missing
%     features now appear as NaN rather than artificial zeros.
%   - Bombcell QC metrics are stored per unit when enabled.

%% 0 - Setup
addpath(genpath("/path/to/DeePhys")) % CHANGE

%% 1 - Load a single recording
% Point to a specific MEArecording.mat file
file_path = "/path/to/sorting/MEArecording.mat"; % CHANGE
load(file_path)

%% 2 - Explore the MEArecording object
% The loaded variable 'obj' is an MEArecording with these properties:
%
%   Metadata        - from the lookup xlsx (Mutation, DIV, ChipID, etc.)
%   RecordingInfo   - sampling rate, duration, electrode coordinates
%   Parameters      - QC thresholds, analysis toggles, connectivity params
%   Spikes          - spike times and unit assignments (post-QC)
%   Units           - array of Unit objects that passed QC
%   Bursts          - burst detection results (raw + pruned)
%   NetworkFeatures - Catch22, Regularity, GraphFeatures, BurstFeatures
%   UnitFeatures    - aggregated single-cell features across all units
%   Connectivity    - CCG and STTC graphs + connection lists

obj

%% 3 - Retrieve data
% Original (pre-QC) spike times
[spike_times, spike_units, hasSpikes] = getSpikeTimes(obj);

% Template waveforms (cached during construction, freed after)
template_matrix = getTemplateMatrix(obj);

% Table of all single-cell features for all good units
unit_table = obj.getUnitFeatures(obj, "all");

% Inspect what features are available
disp(unit_table.Properties.VariableNames')

% New features added in 2026 refactor:
%   CV2InterSpikeInterval  - local ISI variability (Holt et al., 1996)
%   LocalVariation         - LV (Shinomoto et al., 2003)
%   RevisedLocalVariation  - LvR with refractory correction (Shinomoto, 2009)
%   FanoFactor             - var(count)/mean(count) in 100ms bins
%   HalfWidth              - waveform duration at half-maximum trough (ms)

%% 4 - Connectivity inspection
% Look for active neurons
[fr, fr_idx] = maxk(unit_table.FiringRate, 10);

% Check excitatory connections from the CCG method
if ~isempty(obj.Connectivity) && isfield(obj.Connectivity, 'CCG')
    exc_con = obj.Connectivity.CCG.ExcitatoryConnection;

    % Find connections between highly active units
    unit_overlap = ismember(exc_con, fr_idx);
    candidate_con = find(sum(unit_overlap, 2) == 2);

    if ~isempty(candidate_con)
        % Plot the CCG for the first candidate pair
        figure('Color', 'w')
        obj.PlotCCG(exc_con(candidate_con(1), 1), exc_con(candidate_con(1), 2));

        % Visualise the unit templates on the electrode array
        unit_ids = exc_con(candidate_con(1), :);
        wf_cutout = 1:60;
        scale_factor = 1.9; % adjust for electrode spacing
        figure('Color', 'w')
        obj.PlotUnitTemplate(unit_ids, scale_factor, wf_cutout);
    end
end

%% 5 - STTC connectivity
if ~isempty(obj.Connectivity) && isfield(obj.Connectivity, 'STTC')
    empirical_sttc = obj.empiricalSTTCmatrix();

    method = "STTC";
    matrix = "wu"; % weighted undirected
    figure('Color', 'w')
    obj.PlotCommunity(method, matrix);
end

%% 6 - Burst detection
% Default ISIN method (adaptive ISI_N threshold)
obj.detectBursts();

% Check results visually
obj.PlotBurstCheck("hist", cutout=[150, 200], binning=0.1)
obj.PlotBurstCheck("scatter", cutout=[150, 200])

% The inferred parameters are stored in obj.Bursts
disp(obj.Bursts)

% Manual parameter override:
%   obj.detectBursts(N=7, ISI_N=0.01);

% NEW: logISI burst detection (Pasquale et al., 2010)
% Fits a 2-component GMM to log10(ISI) to find a data-driven threshold.
% Best suited for cultures with clear bimodal ISI distributions.
%   obj.detectBursts(method="logISI");

%% 7 - Graph features (new in 2026)
% If connectivity was computed, graph features now include:
if isfield(obj.NetworkFeatures, 'GraphFeatures')
    gf = obj.NetworkFeatures.GraphFeatures;
    disp(gf)
    % New: Modularity (Louvain), SmallWorldness (sigma), RichClub
end

%% 8 - Bombcell QC metrics (if enabled during extraction)
% When params.QC.Bombcell.Enable = true, each unit carries:
%   unit.BombcellMetrics  - table with BC_ prefixed features
%   unit.BombcellType     - 1=good, 2=MUA, 3=non-somatic good, 4=non-somatic MUA
%
% Example:
%   obj.Units(1).BombcellMetrics
%   obj.Units(1).BombcellType

%% 9 - Save changes
% If you made changes (e.g. re-ran burst detection with different params):
obj.Parameters.Save.Overwrite = true;
obj.saveObject();
