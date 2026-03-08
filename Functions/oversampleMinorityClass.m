function [resampled_table, y_resampled, resampled_group, train_idx, test_idx] = ...
        oversampleMinorityClass(input_table, object_group, y_train, train_idx)
% OVERSAMPLEMINORITYCLASS  Balance class sizes by oversampling the minority class.
%
% Duplicates rows from the smaller class(es) in the training set until all classes
% have equal representation. Test rows are left unchanged and appended after the
% resampled training rows.
%
% INPUTS:
%   input_table  - (N x F) feature table
%   object_group - (1 x N) object array (e.g. Unit array) corresponding to table rows
%   y_train      - (1 x N_train) integer class labels for training samples
%   train_idx    - (1 x N) logical index identifying training rows in input_table
%
% OUTPUTS:
%   resampled_table  - Feature table after resampling + test rows appended
%   y_resampled      - Class labels for resampled training rows
%   resampled_group  - Object group after resampling + test rows appended
%   train_idx        - Logical index for training rows in resampled output
%   test_idx         - Logical index for test rows in resampled output

train_rows = input_table(train_idx,:);
test_rows  = input_table(~train_idx,:);
train_grp  = object_group(train_idx);
test_grp   = object_group(~train_idx);

classes    = unique(y_train);
max_count  = max(histcounts(y_train, [classes; classes(end)+1]));

new_rows  = {};
new_y     = [];
new_grp   = {};

for c = classes
    class_idx = y_train == c;
    n_class   = sum(class_idx);
    n_needed  = max_count - n_class;

    class_rows = train_rows(class_idx,:);
    class_grp  = train_grp(class_idx);
    class_y    = y_train(class_idx);

    if n_needed > 0
        extra_idx  = randi(n_class, 1, n_needed);
        class_rows = [class_rows; class_rows(extra_idx,:)]; %#ok<AGROW>
        class_grp  = [class_grp,  class_grp(extra_idx)];   %#ok<AGROW>
        class_y    = [class_y,    class_y(extra_idx)];      %#ok<AGROW>
    end

    new_rows{end+1} = class_rows; %#ok<AGROW>
    new_y           = [new_y, class_y]; %#ok<AGROW>
    new_grp{end+1}  = class_grp; %#ok<AGROW>
end

resampled_train   = vertcat(new_rows{:});
resampled_train_g = horzcat(new_grp{:});
n_train           = height(resampled_train);

resampled_table   = [resampled_train; test_rows];
resampled_group   = [resampled_train_g, test_grp];
y_resampled       = new_y;

train_idx = [true(1, n_train),           false(1, height(test_rows))];
test_idx  = [false(1, n_train),          true(1,  height(test_rows))];
end
