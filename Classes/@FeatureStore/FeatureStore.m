classdef FeatureStore < handle
% FEATURESTORE  Table-centric container for features across a collection of recordings.
%
% The central architectural piece of DeePhys v2. ALL downstream analysis (ML,
% statistics, plotting) works with FeatureStore tables — never with MEArecording
% or RecordingGroup objects.
%
% Three tables:
%   UnitTable       — one row per unit.      Columns: UnitID, RecordingID, metadata, features.
%   RecordingTable  — one row per recording. Columns: RecordingID, metadata, features.
%   MetadataTable   — one row per recording. Columns: RecordingID, all metadata fields.
%
% USAGE:
%   % From new pipeline
%   fs = FeatureStore.fromProcessors([proc1, proc2, ...]);
%
%   % Migration from old objects
%   fs = FeatureStore.fromLegacyRecordings(mearec_array);
%
%   % Save / load
%   fs.save('experiment.mat');
%   fs = FeatureStore.load('experiment.mat');
%
%   % Subsetting
%   fs2 = fs.subsetByMetadata('Mutation', {'WT','HET'});
%
%   % Feature matrices for ML
%   [X, unitIDs] = fs.unitMatrix('all');
%   Y = fs.UnitTable.Mutation;
%   result = Classifier.classify(X, Y, opts);

    properties
        UnitTable       table   % one row per unit:      UnitID, RecordingID, metadata cols, feature cols
        RecordingTable  table   % one row per recording: RecordingID, metadata cols, feature cols
        MetadataTable   table   % one row per recording: RecordingID, all metadata fields
    end

    properties (Access = private)
        MetadataFields  string  % names of metadata columns (for unitMatrix "all" exclusion)
    end

    % =====================================================================
    % Construction
    % =====================================================================
    methods (Static)

        function fs = fromProcessors(proc_array)
        % FROMPROCESSORS  Assemble FeatureStore from RecordingProcessor array.
            arguments
                proc_array RecordingProcessor
            end
            unit_tables      = cell(1, numel(proc_array));
            recording_rows   = cell(1, numel(proc_array));
            metadata_rows    = cell(1, numel(proc_array));

            for i = 1:numel(proc_array)
                p = proc_array(i);
                if isempty(p.SpikeData) || isempty(p.SpikeData.RecordingID)
                    warning('FeatureStore:emptyProcessor', ...
                        'Processor %d has no SpikeData — skipping.', i);
                    continue
                end
                rec_id = p.SpikeData.RecordingID;

                % Unit table: feature table + identity/metadata columns
                if ~isempty(p.UnitFeatureTable) && ~isempty(p.Units)
                    uft = p.UnitFeatureTable;
                    unit_ids = string({p.Units.UnitID})';
                    meta_row = FeatureStore.metadataStructToRow(p.SpikeData.Metadata);
                    meta_rep = repmat(meta_row, height(uft), 1);
                    id_tbl   = table(unit_ids, repmat(string(rec_id), height(uft), 1), ...
                                     'VariableNames', {'UnitID','RecordingID'});
                    unit_tables{i} = [id_tbl, meta_rep, uft];
                end

                % Recording table: network features + identity/metadata
                if ~isempty(p.NetworkFeatureTable)
                    meta_row = FeatureStore.metadataStructToRow(p.SpikeData.Metadata);
                    id_tbl   = table(string(rec_id), 'VariableNames', {'RecordingID'});
                    recording_rows{i} = [id_tbl, meta_row, p.NetworkFeatureTable];
                end

                % Metadata table
                meta_row = FeatureStore.metadataStructToRow(p.SpikeData.Metadata);
                id_tbl   = table(string(rec_id), 'VariableNames', {'RecordingID'});
                metadata_rows{i} = [id_tbl, meta_row];
            end

            fs = FeatureStore.assembleFromCells(unit_tables, recording_rows, metadata_rows);
        end

        function fs = fromLegacyRecordings(rec_array, unit_features, network_features)
        % FROMLEGACYRECORDINGS  Build FeatureStore from an old MEArecording array.
        %
        % Migration path: works alongside existing MEArecording objects without
        % requiring a full reprocess. Features are extracted via FeatureAssembly.
            arguments
                rec_array       MEArecording
                unit_features   string = "all"
                network_features string = "all"
            end

            unit_tables    = cell(1, numel(rec_array));
            recording_rows = cell(1, numel(rec_array));
            metadata_rows  = cell(1, numel(rec_array));

            for i = 1:numel(rec_array)
                obj    = rec_array(i);
                rec_id = string(paramHash(struct('InputPath', char(obj.Metadata.InputPath))));

                % Metadata row
                meta_row = FeatureStore.metadataStructToRow(obj.Metadata);
                id_tbl   = table(rec_id, 'VariableNames', {'RecordingID'});
                metadata_rows{i} = [id_tbl, meta_row];

                % Unit table
                if ~isempty(obj.Units)
                    try
                        uft      = FeatureAssembly.unitFeatures(obj, unit_features);
                        unit_ids = string({obj.Units.StableID})';
                        meta_rep = repmat(meta_row, numel(obj.Units), 1);
                        id_u     = table(unit_ids, repmat(rec_id, numel(obj.Units), 1), ...
                                         'VariableNames', {'UnitID','RecordingID'});
                        unit_tables{i} = [id_u, meta_rep, uft];
                    catch ME
                        warning('FeatureStore:fromLegacy', ...
                            'Could not extract unit features for recording %d: %s', i, ME.message);
                    end
                end

                % Recording table
                try
                    rft    = FeatureAssembly.recordingFeatures(obj, network_features, unit_features, false);
                    id_r   = table(rec_id, 'VariableNames', {'RecordingID'});
                    recording_rows{i} = [id_r, meta_row, rft];
                catch ME
                    warning('FeatureStore:fromLegacy', ...
                        'Could not extract recording features for recording %d: %s', i, ME.message);
                end
            end

            fs = FeatureStore.assembleFromCells(unit_tables, recording_rows, metadata_rows);
        end

        function fs = load(file_path)
        % LOAD  Restore a FeatureStore saved with fs.save().
            arguments
                file_path (1,1) string
            end
            s  = builtin('load', file_path, 'UnitTable', 'RecordingTable', 'MetadataTable', 'MetadataFields');
            fs = FeatureStore();
            fs.UnitTable      = s.UnitTable;
            fs.RecordingTable = s.RecordingTable;
            fs.MetadataTable  = s.MetadataTable;
            if isfield(s, 'MetadataFields')
                fs.MetadataFields = s.MetadataFields;
            else
                fs.MetadataFields = fs.inferMetadataFields();
            end
        end

    end

    % =====================================================================
    % Persistence
    % =====================================================================
    methods

        function save(fs, file_path)
        % SAVE  Write tables to a .mat file. Reload with FeatureStore.load().
            arguments
                fs        FeatureStore
                file_path (1,1) string
            end
            UnitTable      = fs.UnitTable;      %#ok<PROP>
            RecordingTable = fs.RecordingTable;  %#ok<PROP>
            MetadataTable  = fs.MetadataTable;   %#ok<PROP>
            MetadataFields = fs.MetadataFields;  %#ok<PROP>
            builtin('save', file_path, 'UnitTable', 'RecordingTable', 'MetadataTable', 'MetadataFields');
        end

    end

    % =====================================================================
    % Subsetting
    % =====================================================================
    methods

        function fs2 = subset(fs, recording_ids)
        % SUBSET  Return a new FeatureStore containing only the specified recordings.
            arguments
                fs            FeatureStore
                recording_ids string
            end
            mask_u = ismember(fs.UnitTable.RecordingID,      recording_ids);
            mask_r = ismember(fs.RecordingTable.RecordingID, recording_ids);
            mask_m = ismember(fs.MetadataTable.RecordingID,  recording_ids);

            fs2 = FeatureStore();
            fs2.UnitTable      = fs.UnitTable(mask_u, :);
            fs2.RecordingTable = fs.RecordingTable(mask_r, :);
            fs2.MetadataTable  = fs.MetadataTable(mask_m, :);
            fs2.MetadataFields = fs.MetadataFields;
        end

        function fs2 = subsetByMetadata(fs, field_name, values)
        % SUBSETBYMETADATA  Subset recordings by a metadata field value.
        %
        %   fs2 = fs.subsetByMetadata('Mutation', {'WT','HET'})
        %   fs2 = fs.subsetByMetadata('DIV', [14, 21])
            arguments
                fs         FeatureStore
                field_name (1,1) string
                values              % string array, cell, or numeric
            end
            if ~ismember(field_name, string(fs.MetadataTable.Properties.VariableNames))
                error('FeatureStore:subsetByMetadata', ...
                    'Field "%s" not found in MetadataTable.', field_name);
            end

            col = fs.MetadataTable.(field_name);
            if ischar(values) || isstring(values) || iscell(values)
                mask = ismember(string(col), string(values));
            else
                mask = ismember(col, values);
            end
            rec_ids = fs.MetadataTable.RecordingID(mask);
            fs2 = fs.subset(rec_ids);
        end

    end

    % =====================================================================
    % Feature matrix extraction
    % =====================================================================
    methods

        function [X, unit_ids] = unitMatrix(fs, feature_groups, parent_features)
        % UNITMATRIX  Return feature table and unit IDs for ML.
        %
        %   [X, unit_ids] = fs.unitMatrix('all')
        %   [X, unit_ids] = fs.unitMatrix(["ActivityFeatures","ReferenceWaveform"])
        %   [X, unit_ids] = fs.unitMatrix('all', ["ACG"])     % prefer parent ACGs
        %   [X, unit_ids] = fs.unitMatrix('all', "all")        % all parent features
        %
        % X is a table with only feature columns (no identity/metadata cols).
        % Output column names are always child names (ACG1, not Parent_ACG1).
        % unit_ids is a string vector of UnitIDs corresponding to rows.
            arguments
                fs              FeatureStore
                feature_groups  string = "all"
                parent_features string = string.empty
            end
            cols     = fs.selectFeatureCols(fs.UnitTable, feature_groups);
            X        = fs.resolveParentFeatures(fs.UnitTable, cols, parent_features);
            unit_ids = fs.UnitTable.UnitID;
        end

        function [X, rec_ids] = recordingMatrix(fs, feature_groups, parent_features)
        % RECORDINGMATRIX  Return feature table and recording IDs for ML.
        %
        %   [X, rec_ids] = fs.recordingMatrix('all')
        %   [X, rec_ids] = fs.recordingMatrix('all', ["ACG"])
            arguments
                fs              FeatureStore
                feature_groups  string = "all"
                parent_features string = string.empty
            end
            cols    = fs.selectFeatureCols(fs.RecordingTable, feature_groups);
            X       = fs.resolveParentFeatures(fs.RecordingTable, cols, parent_features);
            rec_ids = fs.RecordingTable.RecordingID;
        end

        function [X, culture_ids] = cultureMatrix(fs, identity_keys, grouping_var, grouping_values, normalization, feature_groups)
        % CULTUREMATRIX  Culture-level feature table (one row per culture).
        %
        % Groups recordings by identity_keys (e.g. ["ChipID","PlatingDate"]),
        % selects recordings at each value of grouping_var (e.g. DIV = [7 14 21 28]),
        % and returns a wide table with feature columns suffixed by grouping value.
        %
        % Replaces aggregateCultureFeatureTables + Culture object traversal.
            arguments
                fs              FeatureStore
                identity_keys   string              = ["ChipID","PlatingDate"]
                grouping_var    (1,1) string         = "DIV"
                grouping_values                      = [7, 14, 21, 28]
                normalization   (1,1) string         = ""   % "baseline", "scaled", or ""
                feature_groups  string               = "all"
            end

            % Get recording-level feature columns
            [rec_feat, rec_ids] = fs.recordingMatrix(feature_groups);
            meta = fs.MetadataTable;

            % Build culture ID string per recording
            culture_ids_per_rec = fs.buildCultureIDs(meta, identity_keys);

            unique_cultures = unique(culture_ids_per_rec, 'stable');
            wide_rows = {};
            out_culture_ids = string.empty;

            for c = 1:numel(unique_cultures)
                cid = unique_cultures(c);
                % Recordings belonging to this culture
                cult_mask = culture_ids_per_rec == cid;
                cult_rec_ids = meta.RecordingID(cult_mask);

                % For each desired grouping value, find the matching recording
                row_data = {};
                has_all = true;
                for v = 1:numel(grouping_values)
                    val = grouping_values(v);
                    % Find recording in this culture with grouping_var == val (±1 tolerance)
                    gv_col = meta.(grouping_var)(cult_mask);
                    if isnumeric(val)
                        match_mask = abs(gv_col - val) <= 1;
                    else
                        match_mask = string(gv_col) == string(val);
                    end
                    if ~any(match_mask)
                        has_all = false;
                        break
                    end
                    match_rid = cult_rec_ids(find(match_mask, 1));
                    feat_row  = rec_feat(rec_ids == match_rid, :);
                    if isempty(feat_row)
                        has_all = false;
                        break
                    end
                    row_data{v} = feat_row; %#ok<AGROW>
                end

                if ~has_all
                    continue
                end

                % Normalize
                if normalization == "baseline" && ~isempty(row_data)
                    baseline = row_data{1}.Variables;
                    for v = 1:numel(row_data)
                        row_data{v}.Variables = row_data{v}.Variables ./ (baseline + eps);
                    end
                elseif normalization == "scaled" && ~isempty(row_data)
                    all_vals = vertcat(row_data{:}).Variables;
                    lo = min(all_vals); hi = max(all_vals);
                    for v = 1:numel(row_data)
                        row_data{v}.Variables = (row_data{v}.Variables - lo) ./ (hi - lo + eps);
                    end
                end

                % Rename columns with grouping value suffix and concatenate
                wide_row_parts = {};
                for v = 1:numel(row_data)
                    rd = row_data{v};
                    suffix = "_" + string(grouping_values(v));
                    rd.Properties.VariableNames = strcat(rd.Properties.VariableNames, char(suffix));
                    wide_row_parts{v} = rd; %#ok<AGROW>
                end
                wide_rows{end+1} = [wide_row_parts{:}]; %#ok<AGROW>
                out_culture_ids(end+1) = cid; %#ok<AGROW>
            end

            if isempty(wide_rows)
                X = table();
                culture_ids = string.empty;
            else
                X = vertcat(wide_rows{:});
                culture_ids = out_culture_ids(:);
            end
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Access = private)

        function cols = selectFeatureCols(fs, tbl, feature_groups)
        % Return column names matching the requested feature groups.
            all_cols  = string(tbl.Properties.VariableNames);
            id_cols   = ["UnitID","RecordingID"];
            meta_cols = fs.MetadataFields;
            exclude   = [id_cols, meta_cols];

            if isscalar(feature_groups) && feature_groups == "all"
                % Exclude Parent_* columns — they are accessed only via resolveParentFeatures
                cols = all_cols(~ismember(all_cols, exclude) & ~startsWith(all_cols, "Parent_"));
                return
            end

            selected = false(1, numel(all_cols));
            for g = 1:numel(feature_groups)
                fg = feature_groups(g);
                switch fg
                    case "ActivityFeatures"
                        names = Unit.returnFeatureNames("act");
                        selected = selected | ismember(all_cols, names);
                    case "WaveformFeatures"
                        % WaveformFeatures columns: check against known names
                        wf_names = ["PeakTroughRatio","AUC","RiseSlope","DecaySlope","HalfWidth","Asymmetry","ZeroCrossings"];
                        selected = selected | ismember(all_cols, wf_names);
                    case "RegularityFeatures"
                        names = Unit.returnFeatureNames("reg");
                        selected = selected | ismember(all_cols, names);
                    case "Catch22"
                        try
                            names = string(GetAllFeatureNames());
                            selected = selected | ismember(all_cols, names);
                        catch
                            selected = selected | startsWith(all_cols, "SC_");  % catch22 prefix fallback
                        end
                    case "GraphFeatures"
                        selected = selected | startsWith(all_cols, "Graph");
                    case "BombcellMetrics"
                        names = Unit.returnFeatureNames("bc");
                        selected = selected | ismember(all_cols, names);
                    case {"ReferenceWaveform","Waveform"}
                        selected = selected | (startsWith(all_cols,"Waveform") & ~ismember(all_cols, exclude));
                    case "ACG"
                        selected = selected | (startsWith(all_cols,"ACG") & ~startsWith(all_cols,"FullACG"));
                    case "FullACG"
                        selected = selected | startsWith(all_cols,"FullACG");
                    otherwise
                        % Treat as a column-name prefix
                        selected = selected | startsWith(all_cols, fg);
                end
            end
            cols = all_cols(selected & ~ismember(all_cols, exclude) & ~startsWith(all_cols, "Parent_"));
        end

        function out_tbl = resolveParentFeatures(~, tbl, child_cols, parent_features)
        % RESOLVEPARENTFEATURES  Build output table, sourcing from Parent_* where requested.
        %
        % Output always has child column names (ACG1, not Parent_ACG1) so feature
        % matrices from datasets with and without parent features are compatible.
        % Fallback: if Parent_X doesn't exist, child X is used silently.
        %
        %   out_tbl = fs.resolveParentFeatures(tbl, cols, ["ACG"])
        %   out_tbl = fs.resolveParentFeatures(tbl, cols, "all")
        %   out_tbl = fs.resolveParentFeatures(tbl, cols, string.empty)  % child only
            all_cols = string(tbl.Properties.VariableNames);
            out_tbl  = tbl(:, []);   % empty table with same row count
            for i = 1:numel(child_cols)
                col        = child_cols(i);
                parent_col = "Parent_" + col;
                use_parent = ~isempty(parent_features) ...
                    && ismember(parent_col, all_cols) ...
                    && (any(parent_features == "all") || ...
                        FeatureStore.isParentGroup(col, parent_features));
                if use_parent
                    out_tbl.(col) = tbl.(parent_col);
                else
                    out_tbl.(col) = tbl.(col);
                end
            end
        end

        function mf = inferMetadataFields(fs)
        % Infer metadata field names from MetadataTable (all cols except RecordingID).
            all_cols = string(fs.MetadataTable.Properties.VariableNames);
            mf = all_cols(all_cols ~= "RecordingID");
        end

    end

    methods (Static, Access = private)

        function tf = isParentGroup(col_name, parent_groups)
        % ISPARENTGROUP  True if col_name belongs to one of the requested parent groups.
        %   Maps column naming conventions to feature group strings.
            tf = false;
            for g = 1:numel(parent_groups)
                pg = parent_groups(g);
                switch pg
                    case "ACG"
                        tf = startsWith(col_name, "ACG") & ~startsWith(col_name, "FullACG");
                    case "FullACG"
                        tf = startsWith(col_name, "FullACG");
                    case {"ReferenceWaveform","Waveform","WaveformFeatures"}
                        tf = startsWith(col_name, "Waveform");
                    case "ActivityFeatures"
                        try
                            names = Unit.returnFeatureNames("act");
                            tf = ismember(col_name, names);
                        catch
                            tf = false;
                        end
                    case "RegularityFeatures"
                        try
                            names = Unit.returnFeatureNames("reg");
                            tf = ismember(col_name, names);
                        catch
                            tf = false;
                        end
                    otherwise
                        tf = startsWith(col_name, pg);
                end
                if tf, return, end
            end
        end

    end

    methods (Static)

        function culture_ids = getCultureIDsForUnits(unit_rec_ids, meta_table, identity_keys)
        % GETCULTUREIDSFORUNITS  Map each unit's RecordingID to a culture ID string.
        %
        %   culture_ids = FeatureStore.getCultureIDsForUnits(unitTable.RecordingID, ...
        %                     fs.MetadataTable, ["ChipID","PlatingDate"])
        %
        % Returns a string vector the same length as unit_rec_ids.
            culture_ids_per_rec = FeatureStore.buildCultureIDs(meta_table, identity_keys);
            [~, loc] = ismember(string(unit_rec_ids), string(meta_table.RecordingID));
            culture_ids = repmat("", numel(unit_rec_ids), 1);
            valid = loc > 0;
            culture_ids(valid) = culture_ids_per_rec(loc(valid));
        end

        function culture_ids = buildCultureIDs(meta_table, identity_keys)
        % BUILDCULTUREIDS  Build a string culture ID per recording row.
        %
        %   culture_ids = FeatureStore.buildCultureIDs(fs.MetadataTable, ["ChipID","PlatingDate"])
        %
        % Concatenates the values of identity_keys columns with "_" separator.
        % Returns a (N_recordings x 1) string vector.
            n = height(meta_table);
            parts = cell(n, numel(identity_keys));
            for k = 1:numel(identity_keys)
                key = identity_keys(k);
                if ismember(key, string(meta_table.Properties.VariableNames))
                    col = meta_table.(key);
                    parts(:, k) = cellstr(string(col));
                else
                    parts(:, k) = repmat({'_'}, n, 1);
                end
            end
            culture_ids = string(join(parts, '_'));
        end

    end

    methods (Static, Access = private)

        function fs = assembleFromCells(unit_cells, recording_cells, metadata_cells)
        % Vertically concatenate cell arrays of tables (handling missing rows).
            unit_cells      = unit_cells(~cellfun(@isempty, unit_cells));
            recording_cells = recording_cells(~cellfun(@isempty, recording_cells));
            metadata_cells  = metadata_cells(~cellfun(@isempty, metadata_cells));

            fs = FeatureStore();
            if ~isempty(unit_cells)
                fs.UnitTable = FeatureStore.stackTables(unit_cells);
            else
                fs.UnitTable = table();
            end
            if ~isempty(recording_cells)
                fs.RecordingTable = FeatureStore.stackTables(recording_cells);
            else
                fs.RecordingTable = table();
            end
            if ~isempty(metadata_cells)
                fs.MetadataTable = FeatureStore.stackTables(metadata_cells);
            else
                fs.MetadataTable = table();
            end

            % Infer metadata field names
            if ~isempty(fs.MetadataTable) && ~isempty(fs.MetadataTable.Properties.VariableNames)
                all_cols = string(fs.MetadataTable.Properties.VariableNames);
                fs.MetadataFields = all_cols(all_cols ~= "RecordingID");
            else
                fs.MetadataFields = string.empty;
            end
        end

        function T = stackTables(tbl_cells)
        % Vertically concatenate tables that may have different columns.
        % Missing columns are filled with NaN (numeric) or "" (string/categorical),
        % with type inferred by peeking at a table that already has the column.
            var_cells = cellfun(@(t) string(t.Properties.VariableNames)', tbl_cells, 'un', 0);
            all_vars  = unique(vertcat(var_cells{:}));
            padded    = cellfun(@(t) FeatureStore.padTable(t, all_vars, tbl_cells), tbl_cells, 'un', 0);
            T = vertcat(padded{:});
        end

        function t = padTable(t, all_vars, ref_tables)
        % Add missing columns to t, with fill value inferred from ref_tables.
            existing = string(t.Properties.VariableNames);
            missing  = all_vars(~ismember(all_vars, existing));
            for i = 1:numel(missing)
                col_name = missing(i);
                fill_str = FeatureStore.isStringColumn(col_name, ref_tables);
                if fill_str
                    t.(col_name) = repmat("", height(t), 1);
                else
                    t.(col_name) = NaN(height(t), 1);
                end
            end
        end

        function tf = isStringColumn(col_name, tables)
        % Return true if col_name holds string/char/categorical data in any ref table.
            tf = false;
            for k = 1:numel(tables)
                tbl = tables{k};
                if ismember(col_name, string(tbl.Properties.VariableNames))
                    col = tbl.(col_name);
                    tf  = isstring(col) || ischar(col) || iscell(col) || iscategorical(col);
                    return
                end
            end
            % Fallback: guess from name — ID columns are strings
            tf = endsWith(col_name, "ID") || ismember(col_name, ...
                ["Mutation","Genotype","Label","Type","Name","Date","Group","Condition"]);
        end

        function meta_row = metadataStructToRow(meta_struct)
        % Convert a metadata struct to a single-row table.
            fns = fieldnames(meta_struct);
            if isempty(fns)
                meta_row = table();
                return
            end
            row_data = {};
            var_names = {};
            for i = 1:numel(fns)
                v = meta_struct.(fns{i});
                % Skip fields that are non-scalar or complex types
                if isstruct(v) || iscell(v) || (isnumeric(v) && ~isscalar(v)) || isa(v, 'datetime')
                    continue
                end
                if ischar(v) || isstring(v)
                    row_data{end+1} = string(v);    %#ok<AGROW>
                elseif isnumeric(v) && isscalar(v)
                    row_data{end+1} = v;             %#ok<AGROW>
                elseif islogical(v) && isscalar(v)
                    row_data{end+1} = double(v);     %#ok<AGROW>
                else
                    continue
                end
                var_names{end+1} = fns{i};           %#ok<AGROW>
            end
            if isempty(var_names)
                meta_row = table();
            else
                meta_row = cell2table(row_data, 'VariableNames', var_names);
            end
        end

    end
end
