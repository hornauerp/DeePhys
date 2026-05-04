classdef GroupedCV
% GROUPEDCV  Hierarchy-aware cross-validation that respects data structure.
%
%   Standard cvpartition splits at the row level, so units from the same
%   recording (or recordings from the same culture) can leak into both
%   train and test sets, inflating accuracy estimates.
%
%   GroupedCV splits at a specified group level (recording or culture),
%   then expands the masks to the object level (unit, recording).
%
%   USAGE:
%     gcv = GroupedCV(group_ids, K);
%     for k = 1:gcv.NumFolds
%         [train_idx, test_idx] = gcv.fold(k);
%         X_train = X(train_idx,:); X_test = X(test_idx,:);
%     end
%
%   FACTORY METHOD:
%     gcv = GroupedCV.byGroups(group_ids, K, Y)
%
%   The partition is stratified: each fold has approximately the same
%   proportion of each class label (if provided).

    properties (SetAccess = private)
        NumFolds        (1,1) double    % Number of CV folds
        GroupIDs                        % (N x 1) group assignment per object
        UniqueGroups                    % Unique group values
        GroupPartition  cvpartition     % Partition at the group level
        ClassLabels                     % (N x 1) optional class labels for stratification
    end

    methods

        function gcv = GroupedCV(group_ids, K, class_labels)
        % GROUPEDCV  Constructor.
        %
        % INPUTS:
        %   group_ids    - (N x 1) group assignment for each object (e.g., recording index
        %                  per unit, or culture index per recording)
        %   K            - number of folds (-1 for leave-one-group-out)
        %   class_labels - (optional) (N x 1) class labels for stratified splitting
        %                  at the group level. If provided, the dominant class per group
        %                  is used for stratification via MATLAB's mode():
        %                    - Ties go to the smallest numeric label (mode() tiebreaking).
        %                    - Intra-group class imbalance is flattened to one label per
        %                      group, so minority-class units within a group do not
        %                      influence fold stratification. This is an inherent
        %                      limitation of group-level CV: per-fold class balance is
        %                      not guaranteed when groups contain mixed classes.
            arguments
                group_ids
                K (1,1) double = 5
                class_labels = []
            end

            gcv.GroupIDs = group_ids(:);
            gcv.UniqueGroups = unique(group_ids);
            n_groups = length(gcv.UniqueGroups);

            if K == -1
                K = n_groups;  % leave-one-group-out
            end
            if K > n_groups && K ~= -1
                warning('GroupedCV:foldReduction', ...
                    'Requested K=%d folds but only %d unique groups available — reducing to %d folds.', ...
                    K, n_groups, n_groups);
            end
            gcv.NumFolds = min(K, n_groups);

            % Build group-level labels for stratified partitioning
            if ~isempty(class_labels)
                gcv.ClassLabels = class_labels(:);
                group_class = zeros(n_groups, 1);
                for g = 1:n_groups
                    mask = (gcv.GroupIDs == gcv.UniqueGroups(g));
                    group_class(g) = mode(class_labels(mask));
                end
                gcv.GroupPartition = cvpartition(group_class, 'KFold', gcv.NumFolds);
            else
                gcv.ClassLabels = [];
                gcv.GroupPartition = cvpartition(n_groups, 'KFold', gcv.NumFolds);
            end
        end

        function [train_idx, test_idx] = fold(gcv, k)
        % FOLD  Return logical train/test masks for fold k, expanded to object level.
        %
        % INPUTS:
        %   k - fold number (1 to NumFolds)
        % OUTPUTS:
        %   train_idx - (N x 1) logical mask for training objects
        %   test_idx  - (N x 1) logical mask for test objects
            arguments
                gcv GroupedCV
                k (1,1) double
            end

            train_groups = gcv.UniqueGroups(gcv.GroupPartition.training(k));
            test_groups  = gcv.UniqueGroups(gcv.GroupPartition.test(k));

            train_idx = ismember(gcv.GroupIDs, train_groups);
            test_idx  = ismember(gcv.GroupIDs, test_groups);
        end

        function summary = describe(gcv)
        % DESCRIBE  Return a summary struct of the CV configuration.
            summary.NumFolds = gcv.NumFolds;
            summary.NumGroups = length(gcv.UniqueGroups);
            summary.NumObjects = length(gcv.GroupIDs);
            summary.Stratified = ~isempty(gcv.ClassLabels);
            sizes = arrayfun(@(g) sum(gcv.GroupIDs == g), gcv.UniqueGroups);
            summary.GroupSizes = sizes;
            summary.MeanGroupSize = mean(sizes);
            summary.MinGroupSize = min(sizes);
            summary.MaxGroupSize = max(sizes);
        end

    end

    methods (Static)

        function gcv = byGroups(group_ids, K, Y)
        % BYGROUPS  Generic factory: split by arbitrary group IDs.
        %
        % This is the table-centric entry point used by Classifier.classify and
        % other v2 analysis methods.  group_ids can be any vector that supports
        % unique() (strings, integers, etc.).
        %
        %   gcv = GroupedCV.byGroups(fs.UnitTable.RecordingID, 5, Y)
        %   gcv = GroupedCV.byGroups(culture_ids, -1, Y)  % leave-one-out
        %
        % INPUTS:
        %   group_ids - (N x 1) group assignment per object
        %   K         - folds (-1 = leave-one-group-out)
        %   Y         - (N x 1) optional class labels for stratification
            arguments
                group_ids
                K (1,1) double = 5
                Y = []
            end
            [~, ~, group_idx] = unique(group_ids(:), 'stable');
            gcv = GroupedCV(group_idx, K, Y);
        end

    end
end
