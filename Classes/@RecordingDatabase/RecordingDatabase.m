classdef RecordingDatabase < handle
% RECORDINGDATABASE  SQLite-based registry for tracking MEArecordings.
%
%   Singleton class that maintains a persistent SQLite database of all
%   processed recordings, their analysis status, QC results, and units.
%   Auto-registers recordings during construction, tracks analysis steps,
%   and provides a query interface for building RecordingGroups.
%
%   USAGE:
%     db = RecordingDatabase.instance();    % Get shared instance
%     db.registerRecording(mearec);          % Register a recording
%     db.updateProcessingStatus(mearec, 'generateUnits', 'completed');
%     T = db.queryRecordings('Mutation', 'WT', 'DIV', [14 21]);
%     db.getSummary();
%
%   Requires: Database Toolbox (R2020b+). Degrades gracefully if absent.
%
%   Database location: <project_root>/deephys_registry.db
%   Override via environment variable DEEPHYS_DB_PATH.

    properties (SetAccess = private)
        DBPath       string    % Full path to the SQLite file
        IsConnected  logical = false
    end

    properties (Access = private)
        Connection            % database connection object
    end

    properties (Constant, Hidden)
        SchemaVersion = 1
    end

    methods (Access = private)

        function db = RecordingDatabase(db_path)
        % Private constructor — use RecordingDatabase.instance() instead.
            arguments
                db_path string = ""
            end

            if db_path == ""
                % Check environment variable first
                env_path = getenv('DEEPHYS_DB_PATH');
                if ~isempty(env_path)
                    db.DBPath = string(env_path);
                else
                    % Default: project root (parent of Classes/)
                    class_dir = fileparts(fileparts(mfilename('fullpath')));
                    project_root = fileparts(class_dir);
                    db.DBPath = fullfile(project_root, "deephys_registry.db");
                end
            else
                db.DBPath = db_path;
            end

            db.connect();
        end

    end

    methods

        function connect(db)
        % CONNECT  Open (or reopen) the SQLite connection.
            try
                % Ensure parent directory exists
                [parent_dir, ~, ~] = fileparts(db.DBPath);
                if ~isfolder(parent_dir)
                    mkdir(parent_dir);
                end

                if isfile(db.DBPath)
                    db.Connection = sqlite(db.DBPath);
                else
                    db.Connection = sqlite(db.DBPath, "create");
                end
                db.IsConnected = true;

                % Enable WAL mode for parfor compatibility (best-effort)
                try
                    db.sqlExecute("PRAGMA journal_mode=WAL");
                    db.sqlExecute("PRAGMA busy_timeout=5000");
                catch
                    % Some drivers open implicit transactions; PRAGMAs are optional
                end

                db.ensureSchema();
            catch ME
                warning('RecordingDatabase:connectionFailed', ...
                    'Could not connect to database: %s', ME.message);
                db.IsConnected = false;
            end
        end

        function close(db)
        % CLOSE  Close the database connection.
            if db.IsConnected && ~isempty(db.Connection)
                try
                    close(db.Connection);
                catch
                end
                db.IsConnected = false;
            end
        end

        function delete(db)
        % DELETE  Destructor — close connection on cleanup.
            db.close();
        end

    end

    methods (Access = private)

        function ensureSchema(db)
        % ENSURESCHEMA  Create tables if they do not exist.
            if ~db.IsConnected; return; end

            db.sqlExecute([...
                "CREATE TABLE IF NOT EXISTS schema_version (" ...
                "  version INTEGER NOT NULL, " ...
                "  applied_at TEXT DEFAULT (datetime('now'))" ...
                ")"]);

            db.sqlExecute([...
                "CREATE TABLE IF NOT EXISTS recordings (" ...
                "  id INTEGER PRIMARY KEY AUTOINCREMENT, " ...
                "  input_path TEXT UNIQUE NOT NULL, " ...
                "  chip_id TEXT, " ...
                "  plate_id TEXT, " ...
                "  well_id TEXT, " ...
                "  plating_date TEXT, " ...
                "  recording_date TEXT, " ...
                "  div REAL, " ...
                "  mutation TEXT, " ...
                "  custom_metadata TEXT, " ...  % JSON blob
                "  n_total_templates INTEGER, " ...
                "  n_good_units INTEGER, " ...
                "  sampling_rate REAL, " ...
                "  duration_s REAL, " ...
                "  registered_at TEXT DEFAULT (datetime('now')), " ...
                "  last_updated TEXT DEFAULT (datetime('now'))" ...
                ")"]);

            db.sqlExecute([...
                "CREATE TABLE IF NOT EXISTS processing_status (" ...
                "  id INTEGER PRIMARY KEY AUTOINCREMENT, " ...
                "  recording_id INTEGER NOT NULL, " ...
                "  analysis_name TEXT NOT NULL, " ...
                "  status TEXT NOT NULL DEFAULT 'pending', " ...
                "  parameter_hash TEXT, " ...
                "  started_at TEXT, " ...
                "  completed_at TEXT, " ...
                "  duration_s REAL, " ...
                "  error_message TEXT, " ...
                "  FOREIGN KEY (recording_id) REFERENCES recordings(id), " ...
                "  UNIQUE(recording_id, analysis_name, parameter_hash)" ...
                ")"]);

            db.sqlExecute([...
                "CREATE TABLE IF NOT EXISTS qc_results (" ...
                "  id INTEGER PRIMARY KEY AUTOINCREMENT, " ...
                "  recording_id INTEGER NOT NULL, " ...
                "  n_total_templates INTEGER, " ...
                "  n_good_units INTEGER, " ...
                "  n_ks_filtered INTEGER, " ...
                "  n_bombcell_filtered INTEGER, " ...
                "  ks_label_enabled INTEGER, " ...
                "  bombcell_enabled INTEGER, " ...
                "  qc_params_hash TEXT, " ...
                "  qc_params_json TEXT, " ...
                "  FOREIGN KEY (recording_id) REFERENCES recordings(id)" ...
                ")"]);

            db.sqlExecute([...
                "CREATE TABLE IF NOT EXISTS units (" ...
                "  id INTEGER PRIMARY KEY AUTOINCREMENT, " ...
                "  recording_id INTEGER NOT NULL, " ...
                "  stable_id TEXT, " ...
                "  template_id INTEGER, " ...
                "  reference_electrode INTEGER, " ...
                "  ks_label TEXT, " ...
                "  bombcell_type TEXT, " ...
                "  is_good INTEGER, " ...
                "  FOREIGN KEY (recording_id) REFERENCES recordings(id)" ...
                ")"]);

            % Check if we need to stamp the schema version
            rows = db.sqlFetch("SELECT COUNT(*) AS cnt FROM schema_version");
            if istable(rows)
                cnt = rows.cnt;
            else
                cnt = rows{1};
            end
            if cnt == 0
                db.sqlExecute(sprintf( ...
                    "INSERT INTO schema_version (version) VALUES (%d)", ...
                    db.SchemaVersion));
            end
        end

        function ensureConnected(db)
        % ENSURECONNECTED  Reconnect if needed.
            if ~db.IsConnected || isempty(db.Connection)
                db.connect();
            end
        end

        function sqlExecute(db, sql)
        % SQLEXECUTE  Version-safe SQL execution (execute vs exec).
            if ismethod(db.Connection, 'execute')
                execute(db.Connection, sql);
            else
                exec(db.Connection, sql);
            end
        end

        function rows = sqlFetch(db, sql)
        % SQLFETCH  Version-safe SQL fetch.
            rows = fetch(db.Connection, sql);
        end

    end

    methods

        function rec_id = getRecordingID(db, input_path)
        % GETRECORDINGID  Look up the recordings.id for a given input_path.
        %   Returns [] if not found.
            rec_id = [];
            if ~db.IsConnected; return; end
            try
                rows = db.sqlFetch( ...
                    sprintf("SELECT id FROM recordings WHERE input_path = '%s'", ...
                    strrep(char(input_path), "'", "''")));
                if ~isempty(rows) && height(rows) > 0
                    rec_id = rows.id(1);
                end
            catch
                rec_id = [];
            end
        end

    end

    methods (Static)

        function db = instance(db_path, do_reset)
        % INSTANCE  Get the singleton RecordingDatabase.
        %   db = RecordingDatabase.instance()            — default path
        %   db = RecordingDatabase.instance(custom_path) — override path
        %   RecordingDatabase.instance("", true)         — reset singleton
            arguments
                db_path string = ""
                do_reset logical = false
            end

            persistent shared_db
            if do_reset
                if ~isempty(shared_db) && isvalid(shared_db)
                    shared_db.close();
                end
                shared_db = [];
                db = [];
                return
            end
            if isempty(shared_db) || ~isvalid(shared_db)
                if db_path == ""
                    shared_db = RecordingDatabase();
                else
                    shared_db = RecordingDatabase(db_path);
                end
            end
            db = shared_db;
        end

        function resetInstance()
        % RESETINSTANCE  Clear the singleton (for testing).
            RecordingDatabase.instance("", true);
        end

        function available = isAvailable()
        % ISAVAILABLE  Check whether the Database Toolbox is licensed.
            available = license('test', 'Database_Toolbox') && ...
                        exist('sqlite', 'class') == 8;
        end

    end
end
