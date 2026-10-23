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

    % Re-link units to their parent MEArecording (stripped by Unit.saveobj
    % to break circular references during serialization).
    if ~isempty(obj.Units)
        for u = 1:length(obj.Units)
            obj.Units(u).MEArecording = obj;
        end
    end

    % --- Migration patches (cumulative, keyed by ClassVersion) ---
    % Each block upgrades from version N-1 to N.
    % Legacy objects without ClassVersion are treated as version 0.

    % Version 2: (re)generate every StableID. The ID includes the sorting
    % identity (units of separately sorted days no longer collide); the ID is a deterministic
    % function of stored fields, so stored IDs of the old scheme are simply replaced.
    if ~isempty(obj.Units)
        for u = 1:length(obj.Units)
            obj.Units(u).generateStableID();
        end
    end

    % Loading has no side effects: recordings are registered in the RecordingDatabase only
    % on request (migrateExistingRecordings, RecordingDatabase.registerRecording), so reading
    % a file (e.g. only its metadata) never writes a registry.
end
