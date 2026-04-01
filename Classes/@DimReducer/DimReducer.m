classdef DimReducer
% DIMREDUCER  Static dimensionality reduction methods on plain feature tables.
%
% Zero coupling to MEArecording, Unit, or RecordingGroup.
%
% USAGE:
%   opts.Method            = 'UMAP';   % 'UMAP', 'PCA', 'tSNE'
%   opts.NDims             = 2;
%   opts.NNneighbors       = 15;
%   opts.NormalizationGroups = fs.UnitTable.PlatingDate;
%   opts.ColorFilePath     = '/path/to/umap/colorsByName.properties';
%
%   result = DimReducer.reduce(X, opts);

    methods (Static)

        function result = reduce(X, opts)
        % REDUCE  Dimensionality reduction on a feature table.
        %
        % INPUTS:
        %   X    - (N x F) table of feature values
        %   opts - struct with optional fields:
        %     .Method            - 'UMAP' (default), 'PCA', 'tSNE'
        %     .NDims             - output dimensions (default 2)
        %     .NNeighbors        - UMAP n_neighbors (default 15)
        %     .MinDist           - UMAP min_dist (default 0.1)
        %     .Spread            - UMAP spread (default 1.0)
        %     .NormalizationGroups - (N x 1) string/numeric for per-group z-score
        %     .ColorFilePath     - path to UMAP colorsByName.properties
        %
        % OUTPUTS:
        %   result - DimReductionResult object
            arguments
                X    table
                opts struct = struct()
            end

            opts = DimReducer.defaultOpts(opts);
            norm_groups = DimReducer.resolveGroups(opts.NormalizationGroups, height(X));

            % Normalize
            norm_data = DimReducer.normalizeFeatures(X, norm_groups);
            feat_names = string(X.Properties.VariableNames);

            t_start = tic;
            switch opts.Method
                case 'UMAP'
                    extra_args = {};
                    if ~isempty(opts.ColorFilePath) && isfile(opts.ColorFilePath)
                        extra_args = {'color_file', opts.ColorFilePath};
                    end
                    [reduction, umap, ~, ~] = run_umap(norm_data, ...
                        'n_components',  opts.NDims, ...
                        'n_neighbors',   opts.NNeighbors, ...
                        'min_dist',      opts.MinDist, ...
                        'spread',        opts.Spread, ...
                        'cluster_detail','adaptive', ...
                        'sgd_tasks',     20, ...
                        'verbose',       'none', ...
                        extra_args{:});
                    graph = umap.search_graph;

                case 'PCA'
                    [~, reduction] = pca(norm_data);
                    graph = [];

                case 'tSNE'
                    reduction = tsne(norm_data, 'Algorithm', 'exact');
                    graph = [];

                otherwise
                    error('DimReducer:reduce', 'Unknown method "%s". Use UMAP, PCA, or tSNE.', opts.Method);
            end

            reduction = reduction(:, 1:opts.NDims);

            result = DimReductionResult( ...
                'Reduction',      reduction, ...
                'ObjectGroup',    [], ...
                'Method',         opts.Method, ...
                'UnitFeatures',   feat_names, ...
                'NetworkFeatures', "", ...
                'Parameters',     struct( ...
                    'n_dims',          opts.NDims, ...
                    'n_neighbors',     opts.NNeighbors, ...
                    'normalization_var', opts.NormalizationGroups));

            if opts.Method == "UMAP" && ~isempty(graph)
                result.Graph = graph;
            end

            elapsed = toc(t_start);
            AnalysisLog.instance().add('DimReducer.reduce', opts, ...
                sprintf('%d objects, %s, %dD', size(reduction, 1), opts.Method, opts.NDims), elapsed);
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Static, Access = private)

        function opts = defaultOpts(opts)
            if ~isfield(opts, 'Method'),             opts.Method = 'UMAP';       end
            if ~isfield(opts, 'NDims'),              opts.NDims = 2;              end
            if ~isfield(opts, 'NNeighbors'),         opts.NNeighbors = 15;        end
            if ~isfield(opts, 'MinDist'),            opts.MinDist = 0.1;          end
            if ~isfield(opts, 'Spread'),             opts.Spread = 1.0;           end
            if ~isfield(opts, 'NormalizationGroups'), opts.NormalizationGroups = []; end
            if ~isfield(opts, 'ColorFilePath'),      opts.ColorFilePath = '';     end
        end

        function g = resolveGroups(groups_input, N)
            if isempty(groups_input)
                g = [];
                return
            end
            if ischar(groups_input) || isstring(groups_input)
                [~, ~, g] = unique(string(groups_input(:)), 'stable');
            else
                [~, ~, g] = unique(groups_input(:), 'stable');
            end
            if numel(g) ~= N
                error('DimReducer:resolveGroups', ...
                    'Group vector length (%d) must match number of observations (%d).', numel(g), N);
            end
        end

        function norm_data = normalizeFeatures(X, norm_groups)
            mat = X.Variables;

            % Impute NaN with column median
            for col = 1:size(mat, 2)
                nan_mask = isnan(mat(:, col));
                if any(nan_mask)
                    med = median(mat(~nan_mask, col));
                    if isnan(med), med = 0; end
                    mat(nan_mask, col) = med;
                end
            end

            if ~isempty(norm_groups)
                np = NormalizationPipeline.groupThenGlobal();
                [norm_data, ~] = np.fit_transform(mat, norm_groups);
            else
                np = NormalizationPipeline.globalOnly();
                [norm_data, ~] = np.fit_transform(mat, []);
            end

            % Remove all-NaN columns
            bad = any(isnan(norm_data));
            norm_data(:, bad) = [];
        end

    end
end
