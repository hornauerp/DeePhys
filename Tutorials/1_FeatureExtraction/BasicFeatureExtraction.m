% BasicFeatureExtraction.m
%
% This tutorial shows how to perform feature extraction from a list of
% Kilosort sorting outputs stored in phy-compatible format
% (https://github.com/cortex-lab/phy).
%
% The pipeline reads spike_times.npy, spike_templates.npy, templates.npy,
% channel_positions.npy, and params.py from each sorting directory, then
% runs: QC -> unit generation -> activity/waveform features -> regularity ->
% Catch22 -> burst detection -> connectivity (CCG + STTC) -> graph features.
%
% WHAT CHANGED (2026 refactor):
%   - Metadata is now handled via a per-recording table instead of the old
%     condition-level lookup with PathPartIndex path parsing. Use
%     generateMetadataTable() to convert old-style lookup sheets.
%   - The MetadataTable is pre-read once and passed via metadata struct,
%     so parfor workers don't each re-read the xlsx from disk.
%   - retrieveMetadata() now supports BOTH formats: per-recording
%     (InputPath column) and legacy condition-level (Chip_IDs column).
%   - Burst detection now supports 'logISI' method (Pasquale et al., 2010)
%     in addition to the default 'ISIN' method.
%   - New activity features: CV2, LV, LvR, Fano factor.
%   - New waveform feature: HalfWidth.
%   - New graph features: Modularity, SmallWorldness, rich club re-enabled.
%   - Bombcell QC and Kilosort label filtering are available via parameters.
%   - If the Database Toolbox is available, recordings are automatically
%     registered in the RecordingDatabase during construction.
%
% OLD APPROACH (pre-refactor):
%   metadata.LookupPath = 'cellline_lookup.xlsx';
%   metadata.PathPartIndex = [8, 9, 10]; % indices into path segments
%   % This relied on a specific folder naming convention where path
%   % segments at fixed positions encoded RecordingDate, PlateID, WellID.
%   % The lookup xlsx had one row per condition with a comma-separated
%   % Chip_IDs column. This was fragile and error-prone.

%% 0 - Setup
% Add DeePhys to the path
addpath(genpath("/path/to/DeePhys")) % CHANGE to your DeePhys root

%% 1 - Generate sorting path list
% Specify where your Kilosort outputs live. The path_logic uses wildcards
% to match the subfolder structure.

root_path = "/path/to/your/data"; % CHANGE
path_logic = {'T00*', 'Network', 'w*', 'sorter_output'};

sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n", length(sorting_path_list))

% If your data is stored differently, you can also provide sorting_path_list
% as a plain string array:
%   sorting_path_list = ["/path/to/sort1", "/path/to/sort2"];

%% 2 - Prepare metadata
% RECOMMENDED: per-recording metadata table (new approach).
%
% If you have an existing condition-level lookup sheet (one row per
% genotype/treatment, Chip_IDs column), convert it first:

condition_lookup = "/path/to/cellline_lookup.xlsx"; % CHANGE
PathPartIndex = [6, 7, 9]; % [RecordingDate, PlateID, WellID] indices

metadata_table = generateMetadataTable(sorting_path_list, condition_lookup, ...
    PathPartIndex, 'OutputPath', 'recording_metadata.xlsx');

% This creates recording_metadata.xlsx with one row per sorting path,
% including all columns from the condition lookup plus InputPath,
% RecordingDate, PlateID, WellID, ChipID, DIV.
%
% You can also create this table manually in Excel -- just make sure it has
% an 'InputPath' column matching your sorting paths.

% Point the metadata struct at the per-recording table
metadata = struct();
metadata.LookupPath = "recording_metadata.xlsx";

% ALTERNATIVE: legacy approach (still supported, not recommended)
%   metadata.LookupPath = condition_lookup;
%   metadata.PathPartIndex = PathPartIndex;

%% 3 - Set parameters
% Default parameters work well for most cases. Here are the most commonly
% adjusted ones. See MEArecording.returnDefaultParams() for the full list.

emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();

% QC thresholds
params.QC.Amplitude             = [];         % [min max] uV, empty = skip
params.QC.FiringRate            = [0.01 100]; % [min max] Hz
params.QC.Axon                  = 0.8;        % amplitude ratio for axonal detection
params.QC.Noise                 = 8;          % sign changes in waveform derivative
params.QC.N_Units               = 10;         % minimum good units to proceed

% Optional: Kilosort label pre-filter
% params.QC.KSLabel.Enable      = true;
% params.QC.KSLabel.AcceptedLabels = "good";  % or ["good", "mua"]

% Optional: Bombcell template-based QC
% params.QC.Bombcell.Enable     = true;
% params.QC.Bombcell.Path       = "/path/to/bombcell"; % if not on MATLAB path

% Analyses to run
params.Analyses.SingleCell      = true;
params.Analyses.Regularity      = true;
params.Analyses.Catch22         = true;
params.Analyses.Bursts          = true;        % set false if cultures dont burst
params.Analyses.Connectivity    = ["STTC", "CCG"]; % takes a long time with many units

% Outlier handling
params.Outlier.Method           = [];          % [] = no outlier removal, "median" for robust

% Saving
params.Save.Flag                = true;        % save MEArecording.mat per sorting
params.Save.Overwrite           = true;

%% 4 - Run feature extraction

parallel = true; % uses parfor -- requires Parallel Computing Toolbox

failed_sortings = generate_MEArecordings_from_sorting_list( ...
    sorting_path_list, metadata, params, parallel);

fprintf("Failed: %d/%d sortings\n", length(failed_sortings), length(sorting_path_list))

%% 5 - Debug failed sortings (optional)
% Re-runs failed sortings without parallelisation for proper error messages
if ~isempty(failed_sortings)
    check_failed_analysis(failed_sortings, metadata, params);
end

%% 6 - Database registration (optional)
% If you have the Database Toolbox, recordings are auto-registered during
% construction. To backfill older recordings:

if RecordingDatabase.isAvailable()
    db = RecordingDatabase.instance();
    db.getSummary();

    % Or migrate existing recordings processed before the database existed:
    % n = migrateExistingRecordings(sorting_path_list);
end
