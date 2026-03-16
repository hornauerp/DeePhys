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
%     gcv = GroupedCV(object_group, group_ids, K);
%     for k = 1:gcv.NumFolds
%         [train_idx, test_idx] = gcv.fold(k);
%         X_train = X(train_idx,:); X_test = X(test_idx,:);
%     end
%
%   FACTORY METHODS:
%     gcv = GroupedCV.byRecording(unit_array, K)
%     gcv = GroupedCV.byCulture(unit_array, K)
%     gcv = GroupedCV.byCulture(recording_array, K)
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
        %                  is used for stratification.
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

        function gcv = byRecording(unit_array, K, class_labels)
        % BYRECORDING  Split at the recording level for unit-level analysis.
        %   Units from the same recording are always in the same fold.
            arguments
                unit_array Unit
                K (1,1) double = 5
                class_labels = []
            end
            recordings = [unit_array.MEArecording];
            [~, ~, group_ids] = unique(recordings);
            gcv = GroupedCV(group_ids, K, class_labels);
        end

        function gcv = byCulture(object_array, K, class_labels)
        % BYCULTURE  Split at the culture level.
        %   All recordings/units from the same culture are in the same fold.
        %   Works with both Unit arrays and MEArecording arrays.
            arguments
                object_array
                K (1,1) double = 5
                class_labels = []
            end
            if isa(object_array, 'Unit')
                recordings = [object_array.MEArecording];
            elseif isa(object_array, 'MEArecording')
                recordings = object_array;
            else
                error('GroupedCV:byCulture', 'Input must be Unit array or MEArecording array.');
            end
            metadata = [recordings.Metadata];
            % Culture identity = ChipID + PlatingDate
            culture_keys = strings(length(recordings), 1);
            for i = 1:length(recordings)
                md = metadata(i);
                chip = "unknown";
                plating = "unknown";
                if isfield(md, 'ChipID'),     chip    = string(md.ChipID);     end
                if isfield(md, 'PlatingDate'), plating = string(md.PlatingDate); end
                culture_keys(i) = chip + "_" + plating;
            end
            if isa(object_array, 'Unit')
                % Expand recording-level culture keys to unit level
                unit_culture_keys = strings(length(object_array), 1);
                for u = 1:length(object_array)
                    rec_idx = find(recordings == object_array(u).MEArecording, 1);
                    unit_culture_keys(u) = culture_keys(rec_idx);
                end
                [~, ~, group_ids] = unique(unit_culture_keys);
            else
                [~, ~, group_ids] = unique(culture_keys);
            end
            gcv = GroupedCV(group_ids, K, class_labels);
        end

    end
end
