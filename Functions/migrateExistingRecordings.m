function n_migrated = migrateExistingRecordings(sorting_paths, mat_filename)
% MIGRATEEXISTINGRECORDINGS  Backfill the RecordingDatabase from saved MAT files.
%
%   n = migrateExistingRecordings(sorting_paths)
%   n = migrateExistingRecordings(sorting_paths, 'MEArecording.mat')
%
%   Iterates over a list of Kilosort sorting directories, loads each
%   MEArecording.mat, and registers it in the database. This populates
%   the registry for recordings that were processed before the database
%   was added.
%
%   INPUTS:
%     sorting_paths  - string array of sorting directories
%     mat_filename   - name of the MAT file (default: 'MEArecording.mat')
%
%   RETURNS:
%     n_migrated - number of recordings successfully registered
%
%   EXAMPLE:
%     paths = ["D:/data/sort1", "D:/data/sort2", "D:/data/sort3"];
%     n = migrateExistingRecordings(paths);

    arguments
        sorting_paths string
        mat_filename string = "MEArecording.mat"
    end

    if ~RecordingDatabase.isAvailable()
        error('migrateExistingRecordings:noToolbox', ...
            'Database Toolbox is required.');
    end

    db = RecordingDatabase.instance();
    n_migrated = 0;
    n_total = length(sorting_paths);

    fprintf('Migrating %d recordings to database...\n', n_total);

    for i = 1:n_total
        mat_path = fullfile(sorting_paths(i), mat_filename);

        if exist(mat_path, 'file') ~= 2
            fprintf('  [%d/%d] SKIP (no MAT file): %s\n', i, n_total, sorting_paths(i));
            continue
        end

        try
            data = load(mat_path);
            if isfield(data, 'obj')
                mearec = data.obj;
            else
                % Try first variable in the file
                vars = fieldnames(data);
                mearec = data.(vars{1});
            end

            if ~isa(mearec, 'MEArecording')
                fprintf('  [%d/%d] SKIP (not MEArecording): %s\n', i, n_total, sorting_paths(i));
                continue
            end

            db.registerRecording(mearec);

            if ~isempty(mearec.Units)
                db.registerUnits(mearec);
                db.updateQCResults(mearec);
            end

            % Mark known completed analyses based on what's populated
            if ~isempty(mearec.Units)
                db.updateProcessingStatus(mearec, 'generateUnits', 'completed');
            end
            if ~isempty(mearec.UnitFeatures)
                db.updateProcessingStatus(mearec, 'aggregateSingleCellFeatures', 'completed');
            end
            if isfield(mearec.NetworkFeatures, 'Catch22') && ~isempty(mearec.NetworkFeatures.Catch22)
                db.updateProcessingStatus(mearec, 'runCatch22', 'completed');
            end
            if ~isempty(mearec.Bursts)
                db.updateProcessingStatus(mearec, 'getBurstStatistics', 'completed');
            end
            if ~isempty(mearec.Connectivity) && ~isempty(fieldnames(mearec.Connectivity))
                db.updateProcessingStatus(mearec, 'inferConnectivity', 'completed');
                if ~isempty(mearec.NetworkFeatures) && isfield(mearec.NetworkFeatures, 'GraphFeatures')
                    db.updateProcessingStatus(mearec, 'inferGraphFeatures', 'completed');
                end
            end

            n_migrated = n_migrated + 1;
            fprintf('  [%d/%d] OK: %s (%d units)\n', i, n_total, sorting_paths(i), length(mearec.Units));

        catch ME
            fprintf('  [%d/%d] ERROR: %s — %s\n', i, n_total, sorting_paths(i), ME.message);
        end
    end

    fprintf('Migration complete: %d/%d recordings registered.\n', n_migrated, n_total);
end
