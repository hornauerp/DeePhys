function summary = getSummary(db)
% GETSUMMARY  Print and return a summary of the recording registry.
%
%   summary = db.getSummary()
%
%   Prints to console:
%     - Total recordings, total units
%     - Recordings by mutation
%     - Processing status breakdown
%     - QC pass rates
%
%   Returns a struct with the aggregate data.

    arguments
        db RecordingDatabase
    end

    summary = struct();
    if ~db.IsConnected
        fprintf('RecordingDatabase: not connected.\n');
        return
    end

    try
        db.ensureConnected();

        % Total recordings
        row = db.sqlFetch("SELECT COALESCE(COUNT(*), 0) AS cnt FROM recordings");
        summary.TotalRecordings = row.cnt(1);

        % Total units
        row = db.sqlFetch("SELECT COALESCE(COUNT(*), 0) AS cnt FROM units");
        summary.TotalUnits = row.cnt(1);

        % Total good units per recording (avg)
        row = db.sqlFetch( ...
            "SELECT COALESCE(AVG(n_good_units), 0) AS avg_units FROM recordings WHERE n_good_units > 0");
        summary.AvgGoodUnits = row.avg_units(1);

        % Print header
        fprintf('\n=== RecordingDatabase Summary ===\n');
        fprintf('  Database: %s\n', db.DBPath);
        fprintf('  Recordings: %d\n', summary.TotalRecordings);
        fprintf('  Units: %d\n', summary.TotalUnits);
        fprintf('  Avg good units/recording: %.1f\n', summary.AvgGoodUnits);

        % By mutation (only if recordings exist)
        if summary.TotalRecordings > 0
            mut_table = db.sqlFetch( ...
                "SELECT mutation, COUNT(*) AS cnt FROM recordings GROUP BY mutation ORDER BY cnt DESC");
            summary.ByMutation = mut_table;

            if ~isempty(mut_table) && height(mut_table) > 0
                fprintf('\n  By mutation:\n');
                for i = 1:height(mut_table)
                    mut_name = mut_table.mutation{i};
                    if isempty(mut_name); mut_name = '<unset>'; end
                    fprintf('    %-20s %d\n', mut_name, mut_table.cnt(i));
                end
            end
        end

        % Processing status (only query if table has rows)
        ps_count = db.sqlFetch("SELECT COALESCE(COUNT(*), 0) AS cnt FROM processing_status");
        if ps_count.cnt(1) > 0
            status_table = db.sqlFetch([...
                "SELECT analysis_name, status, COUNT(*) AS cnt " ...
                "FROM processing_status GROUP BY analysis_name, status " ...
                "ORDER BY analysis_name, status"]);
            summary.ProcessingStatus = status_table;

            if ~isempty(status_table) && height(status_table) > 0
                fprintf('\n  Processing status:\n');
                for i = 1:height(status_table)
                    fprintf('    %-25s %-12s %d\n', ...
                        char(status_table.analysis_name(i)), ...
                        char(status_table.status(i)), ...
                        status_table.cnt(i));
                end
            end
        end

        % QC summary (only query if table has rows)
        qc_count = db.sqlFetch("SELECT COALESCE(COUNT(*), 0) AS cnt FROM qc_results");
        if qc_count.cnt(1) > 0
            qc_table = db.sqlFetch([...
                "SELECT COALESCE(AVG(CAST(n_good_units AS REAL) / NULLIF(n_total_templates, 0)), 0) AS avg_pass_rate, " ...
                "COALESCE(AVG(n_good_units), 0) AS avg_good, COALESCE(AVG(n_total_templates), 0) AS avg_total " ...
                "FROM qc_results"]);
            summary.QCSummary = qc_table;

            if qc_table.avg_pass_rate(1) > 0
                fprintf('\n  QC pass rate: %.1f%% (avg %.0f/%.0f units)\n', ...
                    qc_table.avg_pass_rate(1) * 100, ...
                    qc_table.avg_good(1), ...
                    qc_table.avg_total(1));
            end
        end

        fprintf('================================\n\n');

    catch ME
        warning('RecordingDatabase:summaryFailed', ...
            'Could not generate summary: %s', ME.message);
    end
end
