classdef DimReducer
% DIMREDUCER  Static dimensionality reduction methods on plain feature tables.
%
% Zero coupling to MEArecording, Unit, or RecordingGroup.
%
% USAGE:
%   opts.Method            = 'UMAP';    % 'UMAP', 'PCA', 'tSNE'
%   opts.NDims             = 2;
%   opts.NNeighbors        = 'auto';    % 'auto' = round(sqrt(N)), clamped to [5, 200]
%   opts.Perplexity        = 'auto';    % tSNE only; 'auto' = max(5, min(50, round(N/3)))
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
        %     .NNeighbors        - UMAP n_neighbors; 'auto' (default) = round(sqrt(N))
        %     .Perplexity        - tSNE perplexity; 'auto' (default) = max(5,min(50,N/3))
        %     .MinDist           - UMAP min_dist (default 0.1)
        %     .Spread            - UMAP spread (default 1.0)
        %     .NormalizationGroups - (N x 1) string/numeric for per-group z-score
        %     .ColorFilePath     - path to UMAP colorsByName.properties
        %
        % OUTPUTS:
        %   result - DimReductionResult object
        %
        % AUTO PARAMETER RESOLUTION (resolved after data size N is known):
        %   NNeighbors = 'auto'  → round(sqrt(N)), clamped to [5, 200]
        %   Perplexity = 'auto'  → max(5, min(50, round(N/3)))  [tSNE only]
            arguments
                X    table
                opts struct = struct()
            end

            opts = DimReducer.defaultOpts(opts);
            norm_groups = DimReducer.resolveGroups(opts.NormalizationGroups, height(X));

            % Normalize
            norm_data  = DimReducer.normalizeFeatures(X, norm_groups);
            feat_names = string(X.Properties.VariableNames);
            N          = size(norm_data, 1);

            % Resolve 'auto' parameters now that N is known
            n_neighbors = DimReducer.resolveNNeighbors(opts.NNeighbors, N);
            perplexity  = DimReducer.resolvePerplexity(opts.Perplexity, N);

            if ~isempty(opts.Seed)
                rng(opts.Seed, 'twister');
            end

            t_start = tic;
            switch opts.Method
                case 'UMAP'
                    extra_args = {};
                    if ~isempty(opts.ColorFilePath) && isfile(opts.ColorFilePath)
                        extra_args = {'color_file', opts.ColorFilePath};
                    end
                    [reduction, umap, ~, ~] = run_umap(norm_data, ...
                        'n_components',  opts.NDims, ...
                        'n_neighbors',   n_neighbors, ...
                        'min_dist',      opts.MinDist, ...
                        'spread',        opts.Spread, ...
                        'cluster_detail','adaptive', ...
                        'sgd_tasks',     20, ...
                        'verbose',       'none', ...
                        extra_args{:});
                    graph = umap.search_graph;

                case 'PCA'
                    [~, reduction] = pca(norm_data);
                    max_dims = size(reduction, 2);
                    if opts.NDims > max_dims
                        warning('DimReducer:ndimsClamped', ...
                            'NDims=%d exceeds available PCA components (%d) for N=%d observations. Clamping.', ...
                            opts.NDims, max_dims, N);
                        opts.NDims = max_dims;
                    end
                    graph = [];

                case 'tSNE'
                    tsne_args = {'Algorithm', 'exact', 'Perplexity', perplexity};
                    if ~isempty(opts.Seed)
                        tsne_args = [tsne_args, {'Seed', opts.Seed}];
                    end
                    reduction = tsne(norm_data, tsne_args{:});
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
                    'n_dims',           opts.NDims, ...
                    'n_neighbors',      n_neighbors, ...
                    'perplexity',       perplexity, ...
                    'min_dist',         opts.MinDist, ...
                    'spread',           opts.Spread, ...
                    'normalization_var', opts.NormalizationVar, ...
                    'seed',             opts.Seed));

            if opts.Method == "UMAP" && ~isempty(graph)
                result.Graph = graph;
            end

            elapsed = toc(t_start);
            AnalysisLog.instance().add('DimReducer.reduce', opts, ...
                sprintf('%d objects, %s, %dD, n_neighbors=%d', ...
                    N, opts.Method, opts.NDims, n_neighbors), elapsed);
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Static, Access = private)

        function opts = defaultOpts(opts)
            if ~isfield(opts, 'Method'),             opts.Method = 'UMAP';        end
            if ~isfield(opts, 'NDims'),              opts.NDims = 2;               end
            if ~isfield(opts, 'NNeighbors'),         opts.NNeighbors = 'auto';     end  % resolved to round(sqrt(N)) at runtime
            if ~isfield(opts, 'Perplexity'),         opts.Perplexity = 'auto';     end  % tSNE; resolved to max(5,min(50,round(N/3)))
            if ~isfield(opts, 'MinDist'),            opts.MinDist = 0.1;           end
            if ~isfield(opts, 'Spread'),             opts.Spread = 1.0;            end
            if ~isfield(opts, 'NormalizationGroups'), opts.NormalizationGroups = []; end
            if ~isfield(opts, 'ColorFilePath'),      opts.ColorFilePath = '';      end
            if ~isfield(opts, 'Seed'),               opts.Seed = [];               end
            if ~isfield(opts, 'NormalizationVar'),   opts.NormalizationVar = '';   end
        end

        function n = resolveNNeighbors(val, N)
        % RESOLVENEIGHBORS  Resolve NNeighbors: numeric pass-through or auto = round(sqrt(N)).
        %   User range clamped to [5, 200]; N-1 applied last so UMAP constraint
        %   n_neighbors < N always holds, even when N <= 5.
            if ischar(val) || (isstring(val) && val == "auto")
                n = round(sqrt(N));
            else
                n = round(double(val));
            end
            n = max(5, min(n, 200));
            n = min(n, N - 1);  % UMAP hard requirement — applied after range clamp
            if n < 2
                warning('DimReducer:tooFewPoints', ...
                    'N=%d is too small for meaningful UMAP embedding (n_neighbors=%d).', N, n);
            end
        end

        function p = resolvePerplexity(val, N)
        % RESOLVEPERPLEXITY  Resolve Perplexity: numeric pass-through or auto = round(N/3).
        %   Clamped to [5, 50] and to (N-1)/3 so tSNE constraint is satisfied.
            if ischar(val) || (isstring(val) && val == "auto")
                p = round(N / 3);
            else
                p = round(double(val));
            end
            p = max(5, min(p, min(50, floor((N - 1) / 3))));
        end

        function g = resolveGroups(groups_input, N)
            g = MLUtils.resolveGroups(groups_input, N);
        end

        function norm_data = normalizeFeatures(X, norm_groups)
            [norm_data, ~] = MLUtils.normalizeAll(X, norm_groups);
        end

    end
end
