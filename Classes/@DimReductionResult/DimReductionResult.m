classdef DimReductionResult
% DIMREDUCTIONRESULT  Standalone container for dimensionality reduction outputs.
%
%   Stores reduced coordinates, object group, method details, and parameters
%   — decoupled from RecordingGroup state.

    properties
        Reduction           % (N x D) reduced coordinates
        ObjectGroup         % Corresponding objects (Unit/MEArecording/Culture array)
        Method      string  % "UMAP", "PCA", "tSNE"
        UnitFeatures        % Unit feature groups used (string array)
        NetworkFeatures     % Network feature groups used (string array)
        Graph               % UMAP search graph (if applicable)
        ClusterIndex        % Cluster assignments (if computed)
        Parameters          % Struct of all parameters used
    end

    methods
        function r = DimReductionResult(varargin)
            if nargin > 0
                for i = 1:2:length(varargin)
                    r.(varargin{i}) = varargin{i+1};
                end
            end
        end
    end
end
