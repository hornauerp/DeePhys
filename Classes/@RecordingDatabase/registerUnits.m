function registerUnits(db, mearec)
% REGISTERUNITS  Bulk-insert unit records for a recording.
%
%   db.registerUnits(mearec) deletes existing unit rows for this recording
%   and re-inserts all current good units. This makes re-registration
%   idempotent (safe to call after re-processing).
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

        % Delete existing units for this recording (idempotent re-registration)
        execute(db.Connection, sprintf( ...
            "DELETE FROM units WHERE recording_id = %d", rec_id));

        if isempty(mearec.Units); return; end

        % Build and execute INSERT for each unit
        % Using explicit transaction for bulk insert performance
        execute(db.Connection, "BEGIN TRANSACTION");
        try
            for u = 1:length(mearec.Units)
                unit = mearec.Units(u);

                stable_id = '';
                if ~isempty(unit.StableID)
                    stable_id = char(unit.StableID);
                end

                template_id = 0;
                if isprop(unit, 'TemplateID') && ~isempty(unit.TemplateID)
                    template_id = unit.TemplateID;
                end

                ref_elec = 0;
                if isprop(unit, 'ReferenceElectrode') && ~isempty(unit.ReferenceElectrode)
                    ref_elec = unit.ReferenceElectrode;
                end

                ks_label = '';
                if isprop(unit, 'KSLabel') && ~isempty(unit.KSLabel)
                    ks_label = char(unit.KSLabel);
                end

                bc_type = '';
                if isprop(unit, 'BombcellType') && ~isempty(unit.BombcellType)
                    bc_type = char(string(unit.BombcellType));
                end

                esc = @(s) strrep(char(string(s)), "'", "''");

                sql = sprintf([...
                    "INSERT INTO units " ...
                    "(recording_id, stable_id, template_id, reference_electrode, " ...
                    "ks_label, bombcell_type, is_good) " ...
                    "VALUES (%d, '%s', %d, %d, '%s', '%s', 1)"], ...
                    rec_id, esc(stable_id), template_id, ref_elec, ...
                    esc(ks_label), esc(bc_type));

                execute(db.Connection, sql);
            end
            execute(db.Connection, "COMMIT");
        catch ME_inner
            execute(db.Connection, "ROLLBACK");
            rethrow(ME_inner);
        end

    catch ME
        warning('RecordingDatabase:registerUnitsFailed', ...
            'Could not register units: %s', ME.message);
    end
end
