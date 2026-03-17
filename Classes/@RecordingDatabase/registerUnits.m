function registerUnits(db, mearec)
% REGISTERUNITS  Bulk-insert unit records for a recording.
%
%   db.registerUnits(mearec) deletes existing unit rows for this recording
%   and re-inserts all current good units using a single batched INSERT.
%   This makes re-registration idempotent (safe to call after re-processing).
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
        db.sqlExecute(sprintf( ...
            'DELETE FROM units WHERE recording_id = %d', rec_id));

        if isempty(mearec.Units); return; end

        esc = @(s) strrep(char(string(s)), '''', '''''');

        % Build a single batched INSERT with all unit rows
        value_rows = cell(1, length(mearec.Units));
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

            value_rows{u} = sprintf( ...
                '(%d, ''%s'', %d, %d, ''%s'', ''%s'', 1)', ...
                rec_id, esc(stable_id), template_id, ref_elec, ...
                esc(ks_label), esc(bc_type));
        end

        % SQLite supports multi-row INSERT: INSERT INTO t VALUES (...), (...), ...
        % Batch in chunks of 500 to stay within SQL length limits
        chunk_size = 500;
        for start = 1:chunk_size:length(value_rows)
            stop = min(start + chunk_size - 1, length(value_rows));
            sql = ['INSERT INTO units ' ...
                '(recording_id, stable_id, template_id, reference_electrode, ' ...
                'ks_label, bombcell_type, is_good) VALUES ' ...
                strjoin(value_rows(start:stop), ', ')];
            db.sqlExecute(sql);
        end

    catch ME
        warning('RecordingDatabase:registerUnitsFailed', ...
            'Could not register units: %s', ME.message);
    end
end
