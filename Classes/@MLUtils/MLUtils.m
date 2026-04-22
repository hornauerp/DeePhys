classdef MLUtils
% MLUTILS  Shared static utilities for ML pipeline classes.
%
% Centralises helper functions used by Classifier, Regressor, and DimReducer
% to avoid code duplication.

    methods (Static)

        function g = resolveGroups(groups_input, N)
        % RESOLVEGROUPS  Normalise a group vector to a (N x 1) double index, or [].
        %
        % Accepts string arrays, numeric vectors, or categorical arrays.
        % Returns [] when input is empty or unrecognised.
        %
        %   g = MLUtils.resolveGroups(fs.UnitTable.RecordingID, height(X))
        %   g = MLUtils.resolveGroups([], height(X))   % → []
            if isempty(groups_input)
                g = [];
                return
            end
            if ischar(groups_input) || isstring(groups_input)
                [~, ~, g] = unique(string(groups_input(:)), 'stable');
            elseif isnumeric(groups_input) || iscategorical(groups_input)
                [~, ~, g] = unique(groups_input(:), 'stable');
            else
                g = [];
                return
            end
            if numel(g) ~= N
                error('MLUtils:resolveGroups', ...
                    'Group vector length (%d) must match number of observations (%d).', ...
                    numel(g), N);
            end
        end

    end
end
