function T = queryRecordings(db, options)
% QUERYRECORDINGS  Query the recordings table with flexible filters.
%
%   T = db.queryRecordings()                           — all recordings
%   T = db.queryRecordings('Mutation', 'WT')           — by mutation
%   T = db.queryRecordings('DIV', [14 21])             — DIV range
%   T = db.queryRecordings('ChipID', 'plate1_A1')      — by chip
%   T = db.queryRecordings('PlateID', 'plate1')        — by plate
%   T = db.queryRecordings('Status', 'completed')      — fully processed
%   T = db.queryRecordings('Analysis', 'generateUnits') — has specific analysis
%
%   Returns a MATLAB table with all columns from the recordings table.
%   Returns an empty table if the database is not connected or query fails.

    arguments
        db RecordingDatabase
        options.Mutation string = ""
        options.ChipID string = ""
        options.PlateID string = ""
        options.WellID string = ""
        options.DIV double = []       % scalar or [min max] range
        options.RecordingDate string = ""
        options.Status string = ""    % filter by processing_status
        options.Analysis string = ""  % filter by specific analysis name
        options.Limit (1,1) double = 0
    end

    T = table();
    if ~db.IsConnected; return; end

    try
        db.ensureConnected();

        esc = @(s) strrep(char(string(s)), "'", "''");
        conditions = {};
        needs_join = false;

        if options.Mutation ~= ""
            conditions{end+1} = sprintf("r.mutation = '%s'", esc(options.Mutation));
        end
        if options.ChipID ~= ""
            conditions{end+1} = sprintf("r.chip_id = '%s'", esc(options.ChipID));
        end
        if options.PlateID ~= ""
            conditions{end+1} = sprintf("r.plate_id = '%s'", esc(options.PlateID));
        end
        if options.WellID ~= ""
            conditions{end+1} = sprintf("r.well_id = '%s'", esc(options.WellID));
        end
        if options.RecordingDate ~= ""
            conditions{end+1} = sprintf("r.recording_date = '%s'", esc(options.RecordingDate));
        end

        if ~isempty(options.DIV)
            if isscalar(options.DIV)
                conditions{end+1} = sprintf("r.div = %.4g", options.DIV);
            elseif length(options.DIV) == 2
                conditions{end+1} = sprintf("r.div >= %.4g AND r.div <= %.4g", ...
                    options.DIV(1), options.DIV(2));
            end
        end

        % Join-based filters
        if options.Status ~= ""
            needs_join = true;
            conditions{end+1} = sprintf("ps.status = '%s'", esc(options.Status));
        end
        if options.Analysis ~= ""
            needs_join = true;
            conditions{end+1} = sprintf("ps.analysis_name = '%s'", esc(options.Analysis));
        end

        % Build query
        if needs_join
            sql = "SELECT DISTINCT r.* FROM recordings r " + ...
                  "LEFT JOIN processing_status ps ON r.id = ps.recording_id";
        else
            sql = "SELECT * FROM recordings r";
        end

        if ~isempty(conditions)
            sql = sql + " WHERE " + strjoin(string(conditions), " AND ");
        end

        sql = sql + " ORDER BY r.last_updated DESC";

        if options.Limit > 0
            sql = sql + sprintf(" LIMIT %d", options.Limit);
        end

        T = db.sqlFetch( char(sql));

        % Ensure output is a table
        if ~istable(T)
            T = cell2table(T);
        end

    catch ME
        warning('RecordingDatabase:queryFailed', ...
            'Query failed: %s', ME.message);
        T = table();
    end
end
