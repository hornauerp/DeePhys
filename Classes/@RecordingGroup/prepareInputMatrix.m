function [norm_train, norm_test, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx)
            arguments
                rg RecordingGroup
                input_table table
                object_group %Units/recordings/cultures corresponding to the rows of feat_mat
                normalization_var string = "PlatingDate" %Has to correspond to a MEArecording metadata field
                train_idx (1,:) logical = logical(ones(1,size(input_table,1))) %Default to unsupervised; should be logical and length of object_group
                test_idx (1,:) logical = logical(zeros(1,size(input_table,1)))
            end

            input_mat = input_table.Variables;
            feature_names = string(input_table.Properties.VariableNames);
            input_mat(isnan(input_mat)) = 0; %Handle rare case where NaN appears
            if isempty(normalization_var)
                [norm_train, mean_train, sd_train] = normalize(input_mat(train_idx,:));
                nan_idx = any(isnan(norm_train));
                norm_train(:,nan_idx) = [];
                feature_names(nan_idx) = [];
                scale_factor = max(abs(norm_train)); %Make sure we go below the scale of UMAP
                norm_train = norm_train./scale_factor; %Scale data
                if sum(test_idx) > 0

                    norm_test = normalize(input_mat(test_idx,:),'center',mean_train,'scale',sd_train);
                    norm_test(:,nan_idx) = [];
                    norm_test = norm_test./scale_factor;
                else
                    norm_test = [];
                end
            else
                [norm_train, norm_test] = normalizeByGroup(rg, input_mat, object_group, normalization_var, train_idx, test_idx); %Normalize
                nan_idx = any(isnan(norm_train)) | any(isnan(norm_test),1);
                norm_train(:,nan_idx) = [];
                feature_names(nan_idx) = [];
                scale_factor = max(abs(norm_train)); %Compute scale factor before applying
                scale_factor(scale_factor == 0) = 1;
                norm_train = norm_train./scale_factor; %Scale data
                if sum(test_idx) > 0
                    norm_test(:,nan_idx) = [];%Remove peak NaNs
                    norm_test = norm_test./scale_factor; %Scale data using same factor
                else
                    norm_test = [];
                end
            end

end

function [norm_train_data, norm_test_data] = normalizeByGroup(rg, feat_mat, object_group, grouping_var, train_idx, test_idx)
            arguments
                rg RecordingGroup
                feat_mat {isnumeric}
                object_group %Units/recordings/cultures corresponding to the rows of feat_mat
                grouping_var string = "PlatingDate" %Has to correspond to a MEArecording metadata field
                train_idx (1,:) logical = logical(ones(1,size(feat_mat,1))) %Default to unsupervised; should be binary and length of object_group
                test_idx (1,:) logical = logical(zeros(1,size(feat_mat,1)))
            end

            train_idx = logical(train_idx);
            test_idx = logical(test_idx);

            % switch class(object_group)
            %     case 'Unit'
            %         recordings = [object_group.MEArecording];
            %         metadata = [recordings.Metadata];
            %         metadata = [metadata.(grouping_var)];
            %
            %     case 'MEArecording'
            %         metadata = [object_group.Metadata];
            %         metadata = [metadata.(grouping_var)];
            %
            %     case 'cell'
            %         metadata = cellfun(@(x) string(x(1).Metadata.(grouping_var)),object_group);
            % end
            [iG, G] = rg.combineMetadataIndices(object_group, grouping_var);
            g_idx = unique(iG);
            for g = 1:length(g_idx)
                iBatch = (iG == g_idx(g));
                iBatch_train = iBatch & train_idx(:);
                iBatch_test  = iBatch & test_idx(:);
                % Fit normalization on training data within this group only
                [norm_train_batch, batch_mean, batch_sd] = normalize(feat_mat(iBatch_train,:));
                feat_mat(iBatch_train,:) = norm_train_batch;
                % Apply training group stats to test data in the same group
                if any(iBatch_test)
                    feat_mat(iBatch_test,:) = normalize(feat_mat(iBatch_test,:),'center',batch_mean,'scale',batch_sd);
                end
            end
            if length(G) > 1 %Only normalize again if more than 1 group exists
                [norm_train_data, train_mean, train_sd] = normalize(feat_mat(train_idx,:));
                norm_test_data = normalize(feat_mat(test_idx,:),'center',train_mean,'scale',train_sd);
            else
                norm_train_data = feat_mat(train_idx,:);
                norm_test_data = feat_mat(test_idx,:);
            end
            % scale_factor = max(abs(norm_train_data));
            % norm_train_data = norm_train_data./scale_factor;
            % norm_test_data = norm_test_data./scale_factor;
end
