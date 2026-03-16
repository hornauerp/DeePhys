function T = generateMetadataTable(sorting_paths, condition_lookup_path, PathPartIndex, options)
% GENERATEMETADATATABLE  Expand a condition-level lookup into a per-recording metadata table.
%
%   T = generateMetadataTable(sorting_paths, condition_lookup_path, PathPartIndex)
%   T = generateMetadataTable(..., 'OutputPath', 'recording_metadata.xlsx')
%
%   Takes a list of sorting paths and a condition-level lookup xlsx (one row
%   per genotype/treatment with a comma-separated Chip_IDs column), extracts
%   ChipID from each path, joins with the condition table, and produces a
%   per-recording table with one row per sorting path.
%
%   The output table has columns: InputPath, RecordingDate, PlateID, WellID,
%   ChipID, DIV, plus ALL columns from the condition table except Chip_IDs.
%
%   INPUTS:
%     sorting_paths         - string array of sorting output directories
%     condition_lookup_path - path to condition-level xlsx file
%     PathPartIndex         - [date_idx, plate_idx, well_idx] into path segments
%
%   NAME-VALUE OPTIONS:
%     OutputPath - if non-empty, writes the result to this xlsx path
%
%   OUTPUTS:
%     T - table with one row per sorting path
%
%   EXAMPLE:
%     paths = generate_sorting_path_list(root, {'T00*','Network','w*','sorter_output'});
%     T = generateMetadataTable(paths, 'cellline_lookup.xlsx', [6 7 9], ...
%         'OutputPath', 'recording_metadata.xlsx');
%     % Then use: metadata.LookupPath = 'recording_metadata.xlsx'

    arguments
        sorting_paths string
        condition_lookup_path string
        PathPartIndex (1,3) double
        options.OutputPath string = ""
    end

    % Read condition-level table
    cond_table = readtable(condition_lookup_path, "ReadVariableNames", true, "TextType", "string");
    cond_cols = string(cond_table.Properties.VariableNames);

    if ~any(cond_cols == "Chip_IDs")
        error('DeePhys:missingChipIDs', ...
            'Condition lookup table must have a "Chip_IDs" column.');
    end

    % Pre-parse all condition rows' chip IDs for matching
    n_cond = height(cond_table);
    cond_chips = cell(n_cond, 1);
    cond_matched = false(n_cond, 1);
    for r = 1:n_cond
        cond_chips{r} = strtrim(strsplit(string(cond_table.Chip_IDs(r)), ","));
    end

    % Metadata columns to propagate (everything except Chip_IDs)
    meta_cols = cond_cols(cond_cols ~= "Chip_IDs");

    % Build per-recording rows
    n_paths = length(sorting_paths);
    InputPath = strings(n_paths, 1);
    RecordingDate = strings(n_paths, 1);
    PlateID = strings(n_paths, 1);
    WellID = strings(n_paths, 1);
    ChipID = strings(n_paths, 1);
    DIV = nan(n_paths, 1);

    % Initialize dynamic metadata columns
    meta_vals = cell(1, length(meta_cols));
    for c = 1:length(meta_cols)
        col_data = cond_table.(meta_cols(c));
        if isnumeric(col_data)
            meta_vals{c} = nan(n_paths, 1);
        else
            meta_vals{c} = strings(n_paths, 1);
        end
    end

    unmatched_paths = [];

    for i = 1:n_paths
        InputPath(i) = sorting_paths(i);

        % Extract identifiers from path
        try
            ids = MEArecording.extractPathIdentifiers(sorting_paths(i), PathPartIndex);
        catch ME
            warning('DeePhys:pathExtractFailed', ...
                'Could not extract identifiers from path %d: %s — %s', i, sorting_paths(i), ME.message);
            unmatched_paths = [unmatched_paths, i]; %#ok<AGROW>
            continue
        end

        RecordingDate(i) = ids.RecordingDate;
        PlateID(i) = ids.PlateID;
        WellID(i) = ids.WellID;
        ChipID(i) = ids.ChipID;

        % Find matching condition row
        match_row = [];
        for r = 1:n_cond
            if ismember(ids.ChipID, cond_chips{r})
                match_row = r;
                cond_matched(r) = true;
                break
            end
        end

        if isempty(match_row)
            warning('DeePhys:noConditionMatch', ...
                'ChipID "%s" (path %d) not found in condition table.', ids.ChipID, i);
            unmatched_paths = [unmatched_paths, i]; %#ok<AGROW>
            continue
        end

        % Copy condition metadata
        for c = 1:length(meta_cols)
            val = cond_table.(meta_cols(c))(match_row);
            if iscell(val); val = val{1}; end
            meta_vals{c}(i) = val;
        end

        % Compute DIV if PlatingDate is available
        if any(meta_cols == "PlatingDate")
            plating_val = cond_table.PlatingDate(match_row);
            if iscell(plating_val); plating_val = plating_val{1}; end
            rec_dt = parseDateRobust(ids.RecordingDate);
            plat_dt = parseDateRobust(plating_val);
            if ~isnat(rec_dt) && ~isnat(plat_dt)
                DIV(i) = days(rec_dt - plat_dt);
            end
        end
    end

    % Warn about unmatched condition rows
    for r = 1:n_cond
        if ~cond_matched(r)
            warning('DeePhys:unmatchedConditionRow', ...
                'Condition row %d (Chip_IDs: %s) had no matching sorting paths.', ...
                r, cond_table.Chip_IDs(r));
        end
    end

    % Build output table
    T = table(InputPath, RecordingDate, PlateID, WellID, ChipID, DIV);
    for c = 1:length(meta_cols)
        T.(meta_cols(c)) = meta_vals{c};
    end

    % Remove rows that completely failed path extraction
    if ~isempty(unmatched_paths)
        fprintf('Note: %d/%d paths could not be matched to conditions.\n', ...
            length(unmatched_paths), n_paths);
    end

    % Write to file if requested
    if options.OutputPath ~= ""
        writetable(T, options.OutputPath);
        fprintf('Wrote per-recording metadata table (%d rows, %d columns) to: %s\n', ...
            height(T), width(T), options.OutputPath);
    end

    fprintf('Generated metadata table: %d recordings, %d metadata columns.\n', ...
        height(T), width(T));
end
