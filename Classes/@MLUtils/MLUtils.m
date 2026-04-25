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

        function [X_train, X_test, feat_names] = normalizeFeatures(X, train_idx, test_idx, norm_groups)
        % NORMALIZEFEATURES  Impute NaN, normalise, drop bad columns — for train/test splits.
        %
        % Fits the normalization pipeline on training rows, applies to test rows.
        % Used by Classifier and Regressor.
        %
        % INPUTS:
        %   X           - (N x F) feature table (table or numeric matrix)
        %   train_idx   - logical or numeric training index into X
        %   test_idx    - logical or numeric test index into X ([] → X_test returned as [])
        %   norm_groups - (N x 1) numeric group vector for per-group z-score, or []
        %
        % OUTPUTS:
        %   X_train    - (N_train x F') normalised training matrix (numeric)
        %   X_test     - (N_test  x F') normalised test matrix, or [] if test_idx is empty
        %   feat_names - (1 x F') string feature names, bad columns removed
            if istable(X)
                mat        = X.Variables;
                feat_names = string(X.Properties.VariableNames);
            else
                mat        = double(X);
                feat_names = "F" + (1:size(mat, 2));
            end

            % Impute NaN per column using training-set median
            for col = 1:size(mat, 2)
                train_col = mat(train_idx, col);
                med = median(train_col(~isnan(train_col)));
                if isnan(med), med = 0; end
                mat(isnan(mat(:, col)), col) = med;
            end

            % Normalization pipeline fitted on training set
            if ~isempty(norm_groups)
                np      = NormalizationPipeline.groupThenGlobal();
                g_train = norm_groups(train_idx);
            else
                np      = NormalizationPipeline.globalOnly();
                g_train = [];
            end
            [X_train, np] = np.fit_transform(mat(train_idx, :), g_train);

            % Transform test set
            if isempty(test_idx) || (~islogical(test_idx) && isempty(test_idx)) || ...
                    (islogical(test_idx) && ~any(test_idx))
                X_test = [];
            else
                if ~isempty(norm_groups)
                    g_test = norm_groups(test_idx);
                else
                    g_test = [];
                end
                X_test = np.transform(mat(test_idx, :), g_test);
            end

            % Drop zero-variance / all-NaN columns (fitted on training; applied to both)
            bad = any(isnan(X_train));
            if ~isempty(X_test)
                bad = bad | any(isnan(X_test));
                X_test(:, bad) = [];
            end
            X_train(:, bad) = [];
            feat_names(bad) = [];
        end

        function [norm_data, feat_names] = normalizeAll(X, norm_groups)
        % NORMALIZEALL  Impute NaN, normalise all rows, drop bad columns.
        %
        % No train/test split — fits on the full dataset. Used by DimReducer.
        %
        % INPUTS:
        %   X           - (N x F) feature table (table or numeric matrix)
        %   norm_groups - (N x 1) numeric group vector for per-group z-score, or []
        %
        % OUTPUTS:
        %   norm_data  - (N x F') normalised matrix (numeric)
        %   feat_names - (1 x F') string feature names, bad columns removed
            if istable(X)
                mat        = X.Variables;
                feat_names = string(X.Properties.VariableNames);
            else
                mat        = double(X);
                feat_names = "F" + (1:size(mat, 2));
            end

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
            feat_names(bad)   = [];
        end

        function checkClassBalance(Y)
        % CHECKCLASSBALANCE  Warn if majority class is more than 3x larger than minority.
            classes = unique(Y);
            if numel(classes) < 2, return, end
            counts = arrayfun(@(c) sum(Y == c), classes);
            ratio  = max(counts) / min(counts);
            if ratio > 3
                [~, maj_idx] = max(counts);
                [~, min_idx] = min(counts);
                warning('Classifier:imbalanced', ...
                    'Class imbalance detected (%.1f:1 ratio). Majority "%s" (%d) vs minority "%s" (%d). ' ...
                    'Consider opts.Prior = ''uniform'' to weight classes equally.', ...
                    ratio, string(classes(maj_idx)), counts(maj_idx), ...
                    string(classes(min_idx)), counts(min_idx));
            end
        end

    end
end
