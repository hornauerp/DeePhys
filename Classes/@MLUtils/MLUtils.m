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
            labels = MLUtils.groupLabels(groups_input, N);
            if isempty(labels)
                g = [];
            else
                [~, ~, g] = unique(labels, 'stable');
            end
        end

        function [permute_fn, scheme] = labelPermuter(Y, groups)
        % LABELPERMUTER  Label shuffle for permutation tests that respects the CV groups.
        %
        %   [permute_fn, scheme] = MLUtils.labelPermuter(Y, groups)
        %   Y_perm = permute_fn();   % uses the current random number generator
        %
        % scheme "group": the label is constant within each group (e.g. genotype per
        %   culture), so whole groups are relabelled; shuffling single units would break
        %   the group structure and make the null too easy to beat.
        % scheme "within-group": labels vary within groups; they are shuffled within each.
        % scheme "observation": no groups.
            Y = Y(:);
            if isempty(groups)
                scheme = "observation";
                permute_fn = @() Y(randperm(numel(Y)));
                return
            end
            groups = groups(:);
            [~, ~, y_idx] = unique(Y);
            if all(splitapply(@(v) all(v == v(1)), y_idx, groups))
                scheme = "group";
                first_row = splitapply(@(i) i(1), (1:numel(Y))', groups);
                group_y = Y(first_row);
                permute_fn = @() MLUtils.indexBy(group_y(randperm(numel(group_y))), groups);
            else
                scheme = "within-group";
                rows_of = splitapply(@(i) {i}, (1:numel(Y))', groups);
                permute_fn = @() MLUtils.shuffleWithin(Y, rows_of);
            end
        end

        function y = indexBy(v, idx)
            y = v(idx);
        end

        function Y = shuffleWithin(Y, rows_of)
            Y0 = Y;
            for gi = 1:numel(rows_of)
                rows = rows_of{gi};
                Y(rows) = Y0(rows(randperm(numel(rows))));
            end
        end

        function labels = groupLabels(groups_input, N)
        % GROUPLABELS  Group vector as an (N x 1) string column that keeps the group names, or [].
        %
        % Use for normalisation groups: fitted per-group statistics are keyed by these names,
        % so new data from the same plate/chip maps back to them and an unseen batch is
        % recognised as unseen. (Indices from resolveGroups are only valid within one call.)
        % Accepts string, char, cellstr, numeric, logical, categorical and datetime input.
            if isempty(groups_input)
                labels = [];
                return
            end
            if ischar(groups_input)
                groups_input = cellstr(groups_input);
            end
            labels = string(groups_input(:));
            if numel(labels) ~= N
                error('MLUtils:resolveGroups', ...
                    'Group vector length (%d) must match number of observations (%d).', ...
                    numel(labels), N);
            end
        end

        function [X_train, X_test, feat_names, prep] = normalizeFeatures(X, train_idx, test_idx, norm_groups, pipeline_type, covariate_labels)
        % NORMALIZEFEATURES  Impute NaN, normalise, drop bad columns — for train/test splits.
        %
        % Fits the normalization pipeline on training rows, applies to test rows.
        % Used by Classifier and Regressor.
        %
        % INPUTS:
        %   X                - (N x F) feature table (table or numeric matrix)
        %   train_idx        - logical or numeric training index into X
        %   test_idx         - logical or numeric test index into X ([] → X_test returned as [])
        %   norm_groups      - (N x 1) numeric group vector for per-group z-score / ComBat, or []
        %   pipeline_type    - string: 'combat' → combatThenGlobal; else groupThenGlobal/globalOnly
        %   covariate_labels - (N x 1) biological labels for ComBat covariate protection, or []
        %
        % OUTPUTS:
        %   X_train    - (N_train x F') normalised training matrix (numeric)
        %   X_test     - (N_test  x F') normalised test matrix, or [] if test_idx is empty
        %   feat_names - (1 x F') string feature names, bad columns removed
        %   prep       - struct to replay this preprocessing on new data via
        %                MLUtils.applyPreprocessing (stored as MLResult.Preprocessing)
            if nargin < 5 || isempty(pipeline_type)
                pipeline_type = "";
            end
            if nargin < 6
                covariate_labels = [];
            end

            if istable(X)
                mat        = X.Variables;
                feat_names = string(X.Properties.VariableNames);
            else
                mat        = double(X);
                feat_names = "F" + (1:size(mat, 2));
            end

            % Impute NaN per column using training-set median
            all_nan_cols = false(1, size(mat, 2));
            meds = zeros(1, size(mat, 2));
            for col = 1:size(mat, 2)
                train_col = mat(train_idx, col);
                med = median(train_col(~isnan(train_col)));
                if isnan(med)
                    med = 0;
                    all_nan_cols(col) = true;
                end
                meds(col) = med;
                mat(isnan(mat(:, col)), col) = med;
            end
            if any(all_nan_cols)
                warning('MLUtils:allNaNColumn', ...
                    '%d feature column(s) are entirely NaN in training set — imputed with 0. Consider dropping these features.', ...
                    sum(all_nan_cols));
            end

            % Select normalization pipeline
            use_combat = string(pipeline_type) == "combat";
            if use_combat
                np        = NormalizationPipeline.combatThenGlobal();
                g_train   = norm_groups(train_idx);
                cov_train = covariate_labels(train_idx);
            elseif ~isempty(norm_groups)
                np        = NormalizationPipeline.groupThenGlobal();
                g_train   = norm_groups(train_idx);
                cov_train = [];
            else
                np        = NormalizationPipeline.globalOnly();
                g_train   = [];
                cov_train = [];
            end
            [X_train, np] = np.fit_transform(mat(train_idx, :), g_train, CovariateLabels=cov_train);

            % Transform test set
            if isempty(test_idx) || (~islogical(test_idx) && isempty(test_idx)) || ...
                    (islogical(test_idx) && ~any(test_idx))
                X_test = [];
            else
                if use_combat
                    g_test   = norm_groups(test_idx);
                    cov_test = [];   % never use test labels: they are what is being predicted
                elseif ~isempty(norm_groups)
                    g_test   = norm_groups(test_idx);
                    cov_test = [];
                else
                    g_test   = [];
                    cov_test = [];
                end
                X_test = np.transform(mat(test_idx, :), g_test, CovariateLabels=cov_test);
            end

            % Drop zero-variance / all-NaN columns. The mask is decided on training rows only
            % so the feature set does not depend on the test fold and can be replayed on new
            % data; remaining test NaNs (e.g. a constant column within a test group) are set to
            % 0, the training mean after the final global z-score.
            bad = any(isnan(X_train), 1);
            if ~isempty(X_test)
                X_test(:, bad) = [];
                X_test(isnan(X_test)) = 0;
            end
            X_train(:, bad) = [];
            prep = struct('InputFeatures', feat_names, 'Medians', meds, 'Pipeline', np, ...
                'Keep', ~bad, 'UsesGroups', use_combat || ~isempty(norm_groups));
            feat_names(bad) = [];
        end

        function X_out = applyPreprocessing(prep, X, norm_groups)
        % APPLYPREPROCESSING  Replay a fold's fitted preprocessing on new data.
        %
        % prep is the struct returned by normalizeFeatures (stored as result.Preprocessing):
        % training medians impute NaN, the fitted NormalizationPipeline transforms, and the
        % training column mask is applied. X may be a table (columns matched by name) or a
        % numeric matrix with the original column order.
            if istable(X)
                [found, loc] = ismember(prep.InputFeatures, string(X.Properties.VariableNames));
                if ~all(found)
                    error('MLUtils:missingFeatures', 'New data lacks %d training feature(s), e.g. "%s".', ...
                        sum(~found), prep.InputFeatures(find(~found, 1)));
                end
                mat = X{:, loc};
            else
                mat = double(X);
                assert(size(mat, 2) == numel(prep.InputFeatures), 'MLUtils:featureCount', ...
                    'Expected %d columns, got %d.', numel(prep.InputFeatures), size(mat, 2));
            end
            nan_mask = isnan(mat);
            med_mat = repmat(prep.Medians, size(mat, 1), 1);
            mat(nan_mask) = med_mat(nan_mask);
            if prep.UsesGroups
                if isempty(norm_groups)
                    error('MLUtils:normGroupsRequired', ...
                        'The model was trained with group normalisation; pass the normalisation groups of the new data.');
                end
                g = MLUtils.groupLabels(norm_groups, size(mat, 1));
            else
                g = [];
            end
            X_out = prep.Pipeline.transform(mat, g, CovariateLabels=[]);
            X_out = X_out(:, prep.Keep);
            X_out(isnan(X_out)) = 0;
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
            all_nan_cols = false(1, size(mat, 2));
            for col = 1:size(mat, 2)
                nan_mask = isnan(mat(:, col));
                if any(nan_mask)
                    med = median(mat(~nan_mask, col));
                    if isnan(med)
                        med = 0;
                        all_nan_cols(col) = true;
                    end
                    mat(nan_mask, col) = med;
                end
            end
            if any(all_nan_cols)
                warning('MLUtils:allNaNColumn', ...
                    '%d feature column(s) are entirely NaN — imputed with 0. Consider dropping these features.', ...
                    sum(all_nan_cols));
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
                    ['Class imbalance detected (%.1f:1 ratio). Majority "%s" (%d) vs minority "%s" (%d). ' ...
                    'Consider opts.Prior = ''uniform'' to weight classes equally.'], ...
                    ratio, string(classes(maj_idx)), counts(maj_idx), ...
                    string(classes(min_idx)), counts(min_idx));
            end
        end

    end
end
