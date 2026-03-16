function obj = loadobj(s)
% LOADOBJ  Reconstruct MEArecording from saved data.
%
% Handles both legacy files (object saved directly) and future struct-based
% format. Applies any necessary migration patches for schema changes.
%
% MATLAB calls this automatically during load() when the class defines it
% as a Static method.

    if isstruct(s)
        % Future: reconstruct from struct-based format
        obj = MEArecording();
        fnames = fieldnames(s);
        for i = 1:length(fnames)
            try
                obj.(fnames{i}) = s.(fnames{i});
            catch
                warning('MEArecording:loadobj', 'Could not restore field "%s".', fnames{i});
            end
        end
    else
        % Legacy file or current format: s is already an MEArecording
        obj = s;
    end

    % --- Migration patches (cumulative, keyed by ClassVersion) ---
    % Each block upgrades from version N-1 to N.
    % Legacy objects without ClassVersion are treated as version 0.

    % Version 1: backfill StableID on Units that lack it
    if ~isempty(obj.Units)
        for u = 1:length(obj.Units)
            if isempty(obj.Units(u).StableID)
                obj.Units(u).generateStableID();
            end
        end
    end

    % Register loaded recording in database (backfills existing recordings)
    try
        if RecordingDatabase.isAvailable()
            db = RecordingDatabase.instance();
            db.registerRecording(obj);
            if ~isempty(obj.Units)
                db.registerUnits(obj);
                db.updateQCResults(obj);
            end
        end
    catch ME
        warning('MEArecording:loadobjDB', 'Database registration on load: %s', ME.message);
    end
end
