classdef Experiment < handle
% EXPERIMENT  Thin orchestration layer connecting FeatureStore to analysis methods.
%
% Optional — users can call Classifier/Regressor/DimReducer directly with tables.
% Experiment provides a convenient API that assembles the right matrices from
% the FeatureStore and stores results.
%
% USAGE:
%   % From new pipeline
%   procs = RecordingProcessor.loadMany(paths);
%   exp = Experiment.fromProcessors(procs);
%
%   % From legacy objects
%   exp = Experiment.fromLegacyGroup(rg);
%
%   % Analysis (results stored in exp.Results)
%   exp.classify('Unit',      'Mutation',     struct('Algorithm','rf'));
%   exp.classify('Recording', 'Mutation',     struct('KFold',5));
%   exp.reduce('Unit',        struct('Method','UMAP'));
%   exp.regress('Culture',    'Concentration', struct('KFold',5));
%
%   % Access results
%   exp.Results.Classification.Mutation    % ClassificationResult array
%   exp.Results.DimReduction.Unit.UMAP     % DimReductionResult
%
%   % Get feature matrices directly
%   [X, ids] = exp.FeatureStore.unitMatrix('all');

    properties
        FeatureStore    FeatureStore    % All features and metadata as tables
        Processors      RecordingProcessor  % (1 x N) kept for spike time access
        Parameters      struct          % Selection and analysis parameters
        Results         struct          % Stores Classification / Regression / DimReduction results
    end

    % =====================================================================
    % Construction
    % =====================================================================
    methods (Static)

        function exp = fromProcessors(proc_array, parameters)
        % FROMPROCESSORS  Build Experiment from RecordingProcessor array.
            arguments
                proc_array  RecordingProcessor
                parameters  struct = struct()
            end
            exp = Experiment();
            exp.Processors  = proc_array;
            exp.Parameters  = parseStructParameters(Experiment.returnDefaultParams(), parameters);
            exp.FeatureStore = FeatureStore.fromProcessors(proc_array);
            exp.Results = struct('Classification', struct(), ...
                                 'Regression',     struct(), ...
                                 'DimReduction',   struct());
        end

        function exp = fromLegacyGroup(rg, parameters)
        % FROMLEGACYGROUP  Migrate from an old RecordingGroup.
            arguments
                rg          RecordingGroup
                parameters  struct = struct()
            end
            exp = Experiment();
            exp.Processors  = RecordingProcessor.empty();
            exp.Parameters  = parseStructParameters(Experiment.returnDefaultParams(), parameters);
            exp.FeatureStore = FeatureStore.fromLegacyRecordings(rg.Recordings);
            exp.Results = struct('Classification', struct(), ...
                                 'Regression',     struct(), ...
                                 'DimReduction',   struct());
        end

        function exp = fromPaths(processor_paths, parameters)
        % FROMPATHS  Load processors from file paths and build Experiment.
            arguments
                processor_paths     % string array or cell array of paths
                parameters struct = struct()
            end
            procs = RecordingProcessor.loadMany(processor_paths);
            exp = Experiment.fromProcessors(procs, parameters);
        end

    end

    % =====================================================================
    % Analysis convenience methods
    % =====================================================================
    methods

        function result = classify(exp, level, classification_var, opts)
        % CLASSIFY  Classify recordings or units by a metadata variable.
        %
        %   exp.classify('Unit',      'Mutation', opts)
        %   exp.classify('Recording', 'Mutation', opts)
        %   exp.classify('Culture',   'Mutation', opts)  % requires opts.GroupingVar
        %
        % opts fields (all optional):
        %   .Algorithm, .KFold, .NHyper  — see Classifier.classify
        %   .FeatureGroups               — string array for feature selection (default "all")
        %   .CVLevel                     — 'recording' or 'culture' for CV grouping
        %   .NormalizationVar            — metadata field name for per-group normalization
        %   .GroupingVar, .GroupingValues — for Culture-level aggregation
        %   .Normalization               — "baseline" or "scaled" for culture aggregation
            arguments
                exp                 Experiment
                level               (1,1) string = "Recording"
                classification_var  (1,1) string = "Mutation"
                opts                struct = struct()
            end
            opts = Experiment.fillOpts(opts);

            [X, Y, cv_groups, norm_groups] = exp.prepareML(level, classification_var, opts);

            opts.CVGroups          = cv_groups;
            opts.NormalizationGroups = norm_groups;

            result = Classifier.classify(X, Y, opts);
            exp.Results.Classification.(classification_var) = result;
        end

        function result = regress(exp, level, regression_var, opts)
        % REGRESS  Regress a numeric metadata variable.
            arguments
                exp             Experiment
                level           (1,1) string = "Recording"
                regression_var  (1,1) string = "DIV"
                opts            struct = struct()
            end
            opts = Experiment.fillOpts(opts);

            [X, Y, ~, norm_groups] = exp.prepareML(level, regression_var, opts);
            opts.NormalizationGroups = norm_groups;
            if isnumeric(Y)
                Y_num = Y;
            else
                Y_num = str2double(string(Y));
            end

            result = Regressor.regress(X, Y_num, opts);
            exp.Results.Regression.(regression_var) = result;
        end

        function result = reduce(exp, level, opts)
        % REDUCE  Dimensionality reduction at unit, recording, or culture level.
            arguments
                exp     Experiment
                level   (1,1) string = "Unit"
                opts    struct = struct()
            end
            opts = Experiment.fillOpts(opts);
            fg   = opts.FeatureGroups;
            pf   = opts.ParentFeatures;

            switch level
                case "Unit"
                    [X, ~] = exp.FeatureStore.unitMatrix(fg, pf);
                    norm_groups = exp.metadataColumn(exp.FeatureStore.UnitTable, opts.NormalizationVar);
                case "Recording"
                    [X, ~] = exp.FeatureStore.recordingMatrix(fg, pf);
                    norm_groups = exp.metadataColumn(exp.FeatureStore.RecordingTable, opts.NormalizationVar);
                case "Culture"
                    [X, ~] = exp.FeatureStore.cultureMatrix( ...
                        opts.IdentityKeys, opts.GroupingVar, opts.GroupingValues, ...
                        opts.Normalization, fg);
                    norm_groups = [];
                otherwise
                    error('Experiment:reduce', 'Unknown level "%s".', level);
            end

            opts.NormalizationGroups = norm_groups;
            result = DimReducer.reduce(X, opts);
            if ~isfield(exp.Results.DimReduction, level)
                exp.Results.DimReduction.(level) = struct();
            end
            method = opts.Method;
            exp.Results.DimReduction.(level).(method) = result;
        end

    end

    % =====================================================================
    % Static defaults
    % =====================================================================
    methods (Static)

        function params = returnDefaultParams()
            params.Selection.IdentityKeys = ["ChipID","PlatingDate"];
            params.Selection.GroupingVar  = "DIV";
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Access = private)

        function [X, Y, cv_groups, norm_groups] = prepareML(exp, level, label_var, opts)
        % Extract feature matrix X, labels Y, and group vectors for the requested level.
            fg = opts.FeatureGroups;
            pf = opts.ParentFeatures;
            switch level
                case "Unit"
                    [X, ~] = exp.FeatureStore.unitMatrix(fg, pf);
                    tbl    = exp.FeatureStore.UnitTable;
                    if opts.CVLevel == "culture"
                        cv_groups = FeatureStore.getCultureIDsForUnits( ...
                            tbl.RecordingID, exp.FeatureStore.MetadataTable, ...
                            opts.IdentityKeys);
                    else
                        cv_groups = exp.metadataColumn(tbl, 'RecordingID');
                    end
                    norm_groups = exp.metadataColumn(tbl, opts.NormalizationVar);
                    Y = exp.extractLabel(tbl, label_var);
                case "Recording"
                    [X, ~] = exp.FeatureStore.recordingMatrix(fg, pf);
                    tbl    = exp.FeatureStore.RecordingTable;
                    cv_groups   = [];
                    norm_groups = exp.metadataColumn(tbl, opts.NormalizationVar);
                    Y = exp.extractLabel(tbl, label_var);
                case "Culture"
                    [X, culture_ids] = exp.FeatureStore.cultureMatrix( ...
                        opts.IdentityKeys, opts.GroupingVar, opts.GroupingValues, ...
                        opts.Normalization, fg);
                    cv_groups   = [];
                    norm_groups = [];
                    % Y must be provided explicitly (culture-level labels cannot be
                    % inferred automatically — there is no single metadata row per culture).
                    if isfield(opts, 'Y')
                        Y = opts.Y;
                        assert(numel(Y) == height(X), ...
                            'opts.Y length (%d) must match number of cultures (%d).', ...
                            numel(Y), height(X));
                    else
                        % Attempt to look up label_var from MetadataTable for each culture.
                        % This works when label_var is constant within every culture.
                        Y = exp.cultureLabels(culture_ids, label_var, opts.IdentityKeys);
                    end
                otherwise
                    error('Experiment:prepareML', 'Unknown level "%s".', level);
            end
        end

        function col = metadataColumn(~, tbl, field_name)
        % Extract a column from a table, returning [] if the column doesn't exist.
            if isempty(field_name) || strlength(string(field_name)) == 0
                col = [];
                return
            end
            vars = string(tbl.Properties.VariableNames);
            if ismember(string(field_name), vars)
                col = tbl.(field_name);
            else
                col = [];
            end
        end

        function Y = extractLabel(~, tbl, label_var)
        % Extract a label vector from a table column.
            vars = string(tbl.Properties.VariableNames);
            if ~ismember(string(label_var), vars)
                error('Experiment:extractLabel', ...
                    'Column "%s" not found in table. Available: %s', ...
                    label_var, strjoin(vars, ', '));
            end
            Y = tbl.(label_var);
        end

        function Y = cultureLabels(exp, culture_ids, label_var, identity_keys)
        % Look up label_var from MetadataTable for each culture_id.
        % Assumes label_var is constant within a culture (takes first matching row).
            meta = exp.FeatureStore.MetadataTable;
            vars = string(meta.Properties.VariableNames);
            if ~ismember(string(label_var), vars)
                error('Experiment:cultureLabels', ...
                    'Column "%s" not found in MetadataTable. ' + ...
                    'Pass Y explicitly via opts.Y for culture-level classification.', ...
                    label_var);
            end
            % Build culture IDs per recording (same keys as cultureMatrix)
            rec_culture_ids = FeatureStore.buildCultureIDs(meta, identity_keys);
            % For each culture_id, find the first matching recording and read label
            Y = cell(numel(culture_ids), 1);
            for c = 1:numel(culture_ids)
                idx = find(rec_culture_ids == culture_ids(c), 1);
                if isempty(idx)
                    Y{c} = missing;
                else
                    val = meta.(label_var)(idx);
                    if iscell(val)
                        Y{c} = val{1};
                    else
                        Y{c} = val;
                    end
                end
            end
            % Flatten to vector / string array
            if all(cellfun(@isnumeric, Y))
                Y = cell2mat(Y);
            else
                Y = string(Y);
            end
        end

    end

    methods (Static, Access = private)

        function opts = fillOpts(opts)
            if ~isfield(opts, 'FeatureGroups'),   opts.FeatureGroups = "all";          end
            if ~isfield(opts, 'NormalizationVar'), opts.NormalizationVar = "";          end
            if ~isfield(opts, 'CVLevel'),          opts.CVLevel = "recording";         end
            if ~isfield(opts, 'Method'),           opts.Method = "UMAP";               end
            if ~isfield(opts, 'IdentityKeys'),     opts.IdentityKeys = ["ChipID","PlatingDate"]; end
            if ~isfield(opts, 'GroupingVar'),      opts.GroupingVar = "DIV";           end
            if ~isfield(opts, 'GroupingValues'),   opts.GroupingValues = [7,14,21,28]; end
            if ~isfield(opts, 'Normalization'),    opts.Normalization = "";             end
            if ~isfield(opts, 'ParentFeatures'),   opts.ParentFeatures = ["ACG"];       end
        end

    end
end
