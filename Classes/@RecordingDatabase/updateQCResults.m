function updateQCResults(db, mearec)
% UPDATEQCRESULTS  Store QC results for a recording.
%
%   db.updateQCResults(mearec) extracts unit counts and QC parameters
%   from the MEArecording and stores them in the qc_results table.
%   Replaces any existing row for the same recording + QC params hash.
%
%   Non-blocking: errors are caught and emitted as warnings.

    arguments
        db RecordingDatabase
        mearec MEArecording
    end

    if ~db.IsConnected; return; end

    try
        db.ensureConnected();

        rec_id = db.getRecordingID(mearec.Metadata.InputPath);
        if isempty(rec_id); return; end

        % Extract QC parameters
        qc = mearec.Parameters.QC;
        qc_hash = paramHash(qc);
        qc_json = jsonencode(qc);
        esc = @(s) strrep(char(string(s)), "'", "''");

        % Count templates and units
        n_total = 0;
        if isfield(mearec.RecordingInfo, 'N_templates')
            n_total = mearec.RecordingInfo.N_templates;
        end

        n_good = 0;
        if ~isempty(mearec.Units)
            n_good = length(mearec.Units);
        end

        % KS / Bombcell filter counts (if available from QC)
        ks_enabled = qc.KSLabel.Enable;
        bc_enabled = qc.Bombcell.Enable;
        n_ks_filtered = 0;   % not easily recoverable post-hoc; placeholder
        n_bc_filtered = 0;

        % Delete existing row for this recording + hash, then insert
        execute(db.Connection, sprintf( ...
            "DELETE FROM qc_results WHERE recording_id = %d AND qc_params_hash = '%s'", ...
            rec_id, qc_hash));

        sql = sprintf([...
            "INSERT INTO qc_results " ...
            "(recording_id, n_total_templates, n_good_units, " ...
            "n_ks_filtered, n_bombcell_filtered, " ...
            "ks_label_enabled, bombcell_enabled, " ...
            "qc_params_hash, qc_params_json) " ...
            "VALUES (%d, %d, %d, %d, %d, %d, %d, '%s', '%s')"], ...
            rec_id, n_total, n_good, n_ks_filtered, n_bc_filtered, ...
            ks_enabled, bc_enabled, qc_hash, esc(qc_json));

        execute(db.Connection, sql);

    catch ME
        warning('RecordingDatabase:qcUpdateFailed', ...
            'Could not update QC results: %s', ME.message);
    end
end
