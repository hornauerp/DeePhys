function [norm_train, norm_test, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx)
% PREPAREINPUTMATRIX  Build normalized train/test matrices from a feature table.
%
% Internally uses NormalizationPipeline for all normalization steps. The
% pipeline is fitted on training data only, then applied to test data.
%
% INPUTS:
%   rg                - RecordingGroup
%   input_table       - table of features (rows = objects, columns = features)
%   object_group      - Units/recordings/cultures corresponding to rows
%   normalization_var - metadata field for per-group normalization (e.g. "PlatingDate").
%                       Pass "" or empty to skip per-group normalization.
%   train_idx         - logical mask for training rows
%   test_idx          - logical mask for test rows
%
% OUTPUTS:
%   norm_train     - (N_train x F) normalized training matrix
%   norm_test      - (N_test x F) normalized test matrix (or [])
%   feature_names  - (1 x F) string array of surviving feature names

    arguments
        rg RecordingGroup
        input_table table
        object_group
        normalization_var string = "PlatingDate"
        train_idx (1,:) logical = true(1, size(input_table, 1))
        test_idx  (1,:) logical = false(1, size(input_table, 1))
    end

    input_mat = input_table.Variables;
    feature_names = string(input_table.Properties.VariableNames);
    % Impute NaN with per-column median from TRAINING data only (prevents data leakage)
    for col = 1:size(input_mat, 2)
        train_col = input_mat(train_idx, col);
        col_median = median(train_col(~isnan(train_col)));
        if isnan(col_median)
            col_median = 0; % All-NaN column fallback
        end
        nan_mask = isnan(input_mat(:, col));
        if any(nan_mask)
            input_mat(nan_mask, col) = col_median;
        end
    end

    % Build group labels for per-group normalization
    if ~isempty(normalization_var) && strlength(normalization_var) > 0
        [iG, G] = rg.combineMetadataIndices(object_group, normalization_var);
        use_groups = length(G) > 1;
    else
        iG = [];
        use_groups = false;
    end

    % Select pipeline
    if use_groups
        np = NormalizationPipeline.groupThenGlobal();
        group_labels_train = iG(train_idx);
    else
        np = NormalizationPipeline.globalOnly();
        group_labels_train = [];
    end

    % Fit on training data
    [norm_train, np] = np.fit_transform(input_mat(train_idx,:), group_labels_train);

    % Remove NaN features (can arise from zero-variance columns)
    nan_idx = any(isnan(norm_train));
    if sum(test_idx) > 0
        % Also check test-side NaN when group normalization is used
        if use_groups
            group_labels_test = iG(test_idx);
        else
            group_labels_test = [];
        end
        norm_test = np.transform(input_mat(test_idx,:), group_labels_test);
        nan_idx = nan_idx | any(isnan(norm_test), 1);
        norm_test(:, nan_idx) = [];
    else
        norm_test = [];
    end

    norm_train(:, nan_idx) = [];
    feature_names(nan_idx) = [];
end
