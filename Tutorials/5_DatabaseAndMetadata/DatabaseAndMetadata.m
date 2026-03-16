% DatabaseAndMetadata.m
%
% This tutorial covers:
%   1. The new per-recording metadata format
%   2. Converting old condition-level lookup sheets
%   3. Using the RecordingDatabase for tracking and querying
%   4. Building RecordingGroups from database queries
%
% WHAT CHANGED (2026 refactor):
%   - Metadata format: the old condition-level xlsx (one row per genotype
%     with comma-separated Chip_IDs) is still supported but the recommended
%     format is now per-recording (one row per sorting path with InputPath
%     column). This eliminates PathPartIndex path parsing fragility.
%   - RecordingDatabase: SQLite-based singleton that auto-tracks all
%     processed recordings, their analysis status, QC results, and units.
%   - RecordingGroup.fromDatabase: query-based construction without
%     manually listing paths.
%
% OLD APPROACH:
%   % Condition-level lookup with comma-separated Chip_IDs:
%   %   PlatingDate | Mutation | ... | Chip_IDs
%   %   221207      | LRRK2   | ... | M05856_0,M05856_1,...
%   %
%   % Required PathPartIndex = [8, 9, 10] to extract recording date,
%   % plate ID, and well ID from fixed positions in the file path.
%   % Fragile: any folder restructuring broke everything.

%% 0 - Setup
addpath(genpath("/path/to/DeePhys")) % CHANGE

%% =====================================================================
%% PART 1: METADATA MANAGEMENT
%% =====================================================================

%% 1.1 - The new per-recording metadata format
% The recommended metadata table has one row per sorting path:
%
%   InputPath | RecordingDate | PlateID | WellID | ChipID | DIV | PlatingDate | Mutation | ...
%   /path/s1  | 240610        | T002513 | 15     | T..._15| 14  | 221207      | LRRK2    | ...
%   /path/s2  | 240610        | T002513 | 02     | T..._02| 14  | 221207      | WT       | ...
%
% The 'InputPath' column is the key -- it must match the sorting paths.
% All other columns become Metadata fields on the MEArecording object.
% You can add any custom columns (Treatment, Concentration, etc.).

%% 1.2 - Creating the metadata table manually
% The simplest approach: create an xlsx with InputPath + your metadata.
% No path parsing, no naming conventions required.

%% 1.3 - Converting from old condition-level lookup
% If you have an existing lookup sheet (Chip_IDs column), convert it:

sorting_paths = generate_sorting_path_list("/path/to/data", ...
    {'T00*', 'Network', 'w*', 'sorter_output'});

condition_lookup = "/path/to/cellline_lookup.xlsx"; % old format
PathPartIndex = [6, 7, 9]; % [RecordingDate, PlateID, WellID] in path

T = generateMetadataTable(sorting_paths, condition_lookup, PathPartIndex, ...
    'OutputPath', 'recording_metadata.xlsx');

% This reads the old condition-level sheet, extracts ChipID from each path
% using PathPartIndex, matches against Chip_IDs, and writes a per-recording
% table. The output can be reviewed/edited in Excel before use.
%
% Warnings are emitted for:
%   - Paths where ChipID extraction fails
%   - ChipIDs not found in the condition table
%   - Condition rows with no matching sorting paths

disp(T(1:min(5, height(T)), :)) % preview first rows

%% 1.4 - Using metadata in feature extraction
% Point the metadata struct at the per-recording table:
metadata = struct();
metadata.LookupPath = "recording_metadata.xlsx";
% No PathPartIndex needed -- the InputPath column handles matching.

% The legacy approach still works if needed:
%   metadata.LookupPath = condition_lookup;
%   metadata.PathPartIndex = [6, 7, 9];

%% 1.5 - How retrieveMetadata works (both formats)
% When MEArecording reads the lookup table, it auto-detects the format:
%   - If the table has an 'InputPath' column: per-recording format.
%     Finds the row matching the current InputPath, copies all columns.
%   - If the table has a 'Chip_IDs' column: legacy condition-level format.
%     Extracts ChipID from path using PathPartIndex, matches against
%     comma-separated Chip_IDs, copies condition columns.
%
% Both formats support arbitrary additional columns -- everything except
% InputPath/Chip_IDs becomes a Metadata field on the MEArecording.

%% =====================================================================
%% PART 2: RECORDING DATABASE
%% =====================================================================

%% 2.1 - Overview
% RecordingDatabase is a SQLite-based singleton that automatically tracks:
%   - All processed recordings (metadata + recording info)
%   - Per-analysis processing status (running/completed/failed + timing)
%   - QC results (unit counts, filter settings)
%   - Individual units (StableID, TemplateID, KSLabel, BombcellType)
%
% Requires: Database Toolbox (R2020b+). Degrades gracefully if absent.
% Database location: <DeePhys_root>/deephys_registry.db

%% 2.2 - Automatic registration
% Recordings are auto-registered during MEArecording construction:
%   1. After parseRecordingInfo: registers metadata + recording info
%   2. During performAnalyses: tracks each step as running/completed/failed
%   3. After performAnalyses: registers units and QC results
%
% Also auto-registers when loading saved MEArecording.mat files (via loadobj).
% All database calls are non-blocking (try-catch) -- pipeline never crashes.

%% 2.3 - Check if available
if RecordingDatabase.isAvailable()
    disp('Database Toolbox available -- RecordingDatabase is active')
else
    disp('Database Toolbox not found -- RecordingDatabase is disabled')
    disp('The pipeline works normally, just without persistent tracking.')
    return % skip remaining database sections
end

%% 2.4 - Access the singleton
db = RecordingDatabase.instance();

% Custom database path (optional):
%   db = RecordingDatabase.instance("/path/to/custom_registry.db");
% Or via environment variable:
%   setenv('DEEPHYS_DB_PATH', '/path/to/custom_registry.db');

%% 2.5 - Summary
db.getSummary();
% Prints: total recordings, total units, avg units/recording,
% breakdown by mutation, processing status, QC pass rates.

%% 2.6 - Query recordings
% Find all WT recordings:
T_wt = db.queryRecordings('Mutation', 'WT');

% Find recordings in a DIV range:
T_div = db.queryRecordings('DIV', [14 21]);

% Combine filters:
T_filtered = db.queryRecordings('Mutation', 'LRRK2', 'DIV', [14 21]);

% Filter by processing status:
T_done = db.queryRecordings('Status', 'completed');

% Filter by specific analysis:
T_ccg = db.queryRecordings('Analysis', 'inferConnectivity');

disp(T_wt)

%% 2.7 - Build RecordingGroup from database
% Instead of manually listing paths, query the database:
rg = RecordingGroup.fromDatabase(struct('Mutation', 'WT', 'DIV', [14 21]));

% This queries the database, loads MEArecording.mat files, and constructs
% the RecordingGroup. Equivalent to the manual approach but without
% maintaining path lists.

%% 2.8 - Migrate existing recordings
% If you have recordings processed before the database was added:
existing_paths = generate_sorting_path_list("/path/to/old/data", ...
    {'T00*', 'Network', 'w*', 'sorter_output'});

n = migrateExistingRecordings(existing_paths);
fprintf('Migrated %d recordings\n', n)

% This loads each MEArecording.mat and registers it in the database,
% including analysis status inferred from populated fields.

%% 2.9 - Database internals
% Tables: recordings, processing_status, qc_results, units, schema_version
% Journal mode: WAL (safe for parfor -- each worker gets its own connection)
% Busy timeout: 5000 ms (handles lock contention)
% Schema version: tracked for future migrations
%
% The database file is excluded from git via .gitignore:
%   deephys_registry.db
%   deephys_registry.db-wal
%   deephys_registry.db-shm
