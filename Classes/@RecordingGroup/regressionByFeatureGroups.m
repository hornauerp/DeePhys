function result = regressionByFeatureGroups(rg, level, regression_var, stratification_var, stratification_values, pooling_vals, grouping_var, grouping_values, ...
                network_features, unit_features, useClustered, normalization_var, normalization, tolerance, N_hyper, K_fold)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                regression_var string = "Concentration" %Metadata variable that will be regressed for
                stratification_var = "Mutation" %Specify the variable by which to split training and test dataset (e.g. train on wildtype, test on mutation)
                stratification_values = [] %Corresponding to stratification_var, if a value is specified then this will be used as training data
                pooling_vals cell = {}
                grouping_var = "DIV"
                grouping_values = nan
                network_features string = "all"
                unit_features string = "all"
                useClustered logical = false
                normalization_var string = "PlatingDate"
                normalization string = []
                tolerance double = 1
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
            end

            if any(arrayfun(@(x) ~isfield(rg.Recordings(x).Metadata,regression_var),1:length(rg.Recordings)))
                error('Regression data missing')
            end

            [~, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);

            if isempty(stratification_values)
                stratification_var = [stratification_var regression_var];
                [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);
                if K_fold == -1 %set kfold to LOOCV
                    K_fold = length(rec_group_idx);
                end
                cv = cvpartition(rec_group_idx,'KFold',K_fold);
                t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');

                for k = 1:K_fold
                    if level == "Unit"
                        error('Not yet implemented')
%                         train_table = object_group(cv.training(k)).getUnitFeatures(unit_features);
%                         test_table = object_group(cv.test(k)).getUnitFeatures(unit_features);
%                         train_idx = logical([ones(1, size(train_table,1)) zeros(1, size(test_table,1))]);
%                         test_idx = ~train_idx;
%                         input_group = [object_group(cv.training(k)).Units, object_group(cv.test(k)).Units];
%                         input_table = [train_table;test_table];
%                         unit_recordings = [object_group.MEArecording];
%                         metadata = [unit_recordings.Metadata];
%                         true_value = [metadata.(regression_var)];
%                         Y_train = true_value(train_idx);
%                         Y_test = true_value(test_idx);
                    else
                        %                         input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                        [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
                        train_idx = cv.training(k);
                        test_idx = cv.test(k);
                        [reg_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, regression_var, pooling_vals);
                        true_value = group_labels_comb(reg_group_idx);
                        if isstring(true_value)
                           true_value = str2double(true_value);
                        end
                        input_group = object_group;
                        Y_train = true_value(train_idx);
                        Y_test = true_value(test_idx);
                    end

                    [X_train, X_test] = prepareInputMatrix(rg, input_table, input_group, normalization_var, train_idx, test_idx);

                    clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t,'Options',statset('UseParallel',true));
                    Y_pred = predict(clf, X_test);
                    %                     predImp = clf.oobPermutedPredictorImportance;
                    predImp = [];

                    result(k).Mdl = clf;
                    result(k).Y_pred = Y_pred;
                    result(k).Y_test = Y_test;
                    result(k).mse_train = resubLoss(clf);
                    result(k).objects = object_group(test_idx);
                    result(k).predImp = predImp;
                    result(k).GroupLabels = group_labels_comb;
                end
            else
                error('Not yet implemented')
                train_group = find(group_labels_comb == stratification_values);

                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [object_group.Units];
                    [unit_group_idx, ~] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);
                    train_idx = unit_group_idx == train_group;
                    test_idx = unit_group_idx ~= train_group;
                    unit_recordings = [object_group.MEArecording];
                    metadata = [unit_recordings.Metadata];
                    true_value = [metadata.(regression_var)];

                elseif level == "Recording"
                    input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                    train_idx = rec_group_idx == train_group;
                    test_idx = rec_group_idx ~= train_group;
                    metadata = [object_group.Metadata];
                    true_value = [metadata.(regression_var)];

                end
                [X_train, X_test] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);
                Y_train = true_value(train_idx);
                Y_test = true_value(test_idx);
                t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');
                clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t,'Options',statset('UseParallel',true));
                Y_pred = predict(clf,X_test);
                %                     predImp = clf.oobPermutedPredictorImportance;
                predImp = [];

                result(k).Mdl = clf;
                result(k).Y_pred = Y_pred;
                result(k).Y_test = Y_test;
                result(k).mse_train = resubLoss(clf);
                result(k).objects = object_group(test_idx);
                result(k).predImp = predImp;
                result(k).GroupLabels = group_labels_comb;

            end
            rg.Regression.(regression_var) = result;
end
