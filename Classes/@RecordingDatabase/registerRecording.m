function registerRecording(db, mearec)
% REGISTERRECORDING  Insert or update a recording in the registry.
%
%   db.registerRecording(mearec) extracts metadata and recording info from
%   the MEArecording object and upserts a row in the recordings table.
%   Known metadata fields map to named columns; any extra fields are stored
%   as a JSON blob in custom_metadata.
%
%   Non-blocking: errors are caught and emitted as warnings.

    arguments
        db RecordingDatabase
        mearec MEArecording
    end

    if ~db.IsConnected; return; end

    try
        db.ensureConnected();

        input_path = string(mearec.Metadata.InputPath);
        md = mearec.Metadata;
        ri = mearec.RecordingInfo;

        % Extract known fields with safe defaults
        chip_id   = getFieldSafe(md, 'ChipID', '');
        plate_id  = getFieldSafe(md, 'PlateID', '');
        well_id   = getFieldSafe(md, 'WellID', '');
        plating_date    = getFieldSafe(md, 'PlatingDate', '');
        recording_date  = getFieldSafe(md, 'RecordingDate', '');
        div       = getFieldSafe(md, 'DIV', NaN);
        mutation  = getFieldSafe(md, 'Mutation', '');

        % Recording info
        sampling_rate    = getFieldSafe(ri, 'SamplingRate', NaN);
        duration_s       = getFieldSafe(ri, 'Duration', NaN);
        n_total_templates = getFieldSafe(ri, 'N_templates', 0);

        % Good units count
        if ~isempty(mearec.Units)
            n_good_units = length(mearec.Units);
        else
            n_good_units = 0;
        end

        % Collect any extra metadata fields into a JSON blob
        known_fields = ["InputPath", "ChipID", "PlateID", "WellID", ...
            "PlatingDate", "RecordingDate", "DIV", "Mutation", ...
            "PathPartIndex", "LookupPath"];
        md_fields = string(fieldnames(md));
        extra_fields = setdiff(md_fields, known_fields);
        custom = struct();
        for i = 1:length(extra_fields)
            custom.(extra_fields(i)) = md.(extra_fields(i));
        end
        if isempty(fieldnames(custom))
            custom_json = '{}';
        else
            custom_json = jsonencode(custom);
        end

        % Escape single quotes for SQL
        esc = @(s) strrep(char(string(s)), "'", "''");

        sql = sprintf([...
            "INSERT INTO recordings " ...
            "(input_path, chip_id, plate_id, well_id, plating_date, " ...
            "recording_date, div, mutation, custom_metadata, " ...
            "n_total_templates, n_good_units, sampling_rate, duration_s) " ...
            "VALUES ('%s', '%s', '%s', '%s', '%s', '%s', %s, '%s', '%s', %d, %d, %s, %s) " ...
            "ON CONFLICT(input_path) DO UPDATE SET " ...
            "chip_id=excluded.chip_id, plate_id=excluded.plate_id, " ...
            "well_id=excluded.well_id, plating_date=excluded.plating_date, " ...
            "recording_date=excluded.recording_date, div=excluded.div, " ...
            "mutation=excluded.mutation, custom_metadata=excluded.custom_metadata, " ...
            "n_total_templates=excluded.n_total_templates, " ...
            "n_good_units=excluded.n_good_units, " ...
            "sampling_rate=excluded.sampling_rate, " ...
            "duration_s=excluded.duration_s, " ...
            "last_updated=datetime('now')"], ...
            esc(input_path), esc(chip_id), esc(plate_id), esc(well_id), ...
            esc(plating_date), esc(recording_date), ...
            nanStr(div), esc(mutation), esc(custom_json), ...
            n_total_templates, n_good_units, ...
            nanStr(sampling_rate), nanStr(duration_s));

        execute(db.Connection, sql);

    catch ME
        warning('RecordingDatabase:registerFailed', ...
            'Could not register recording: %s', ME.message);
    end
end

function val = getFieldSafe(s, field, default)
% GETFIELDSAFE  Return struct field value or default if missing/empty.
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = default;
    end
end

function s = nanStr(val)
% NANSTR  Convert numeric to SQL-safe string (NULL for NaN).
    if isnan(val)
        s = 'NULL';
    else
        s = sprintf('%.6g', val);
    end
end
