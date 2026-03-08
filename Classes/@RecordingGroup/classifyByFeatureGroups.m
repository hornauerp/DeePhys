function result = classifyByFeatureGroups(rg, level, alg, classification_var, pooling_vals, network_features, unit_features, feature_names, useClustered, ...
                                    grouping_var, grouping_values, normalization, normalization_var, N_hyper, K_fold, tolerance)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                alg string = "rf" %rf, svm,
                classification_var string = "Mutation" %Metadata field that determines y_test/y_train
                pooling_vals cell = {} %Values corresponding to classification_var that is being classified for(e.g. "LNA")
                network_features string = "all"
                unit_features string = "all"
                feature_names string = []
                useClustered logical = false
                grouping_var string = [] %Metadata field that groups recordings to cultures
                grouping_values = nan %Selected values corresponding to grouping_var
                normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
                normalization_var string = []%"PlatingDate" % normalization by each value of normalization_var
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = -1 % number of K-fold CV
                tolerance double = 1
            end

            if isempty(grouping_var) %Get group indices for stratification
                object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            end
            [group_idx, group_labels_comb] = rg.combineMetadataIndices(object_group, classification_var, pooling_vals);


            if K_fold == -1 %set kfold to LOOCV
                K_fold = length(group_idx);
            end
            cv = cvpartition(group_idx,'KFold',K_fold);

            if isempty(grouping_var)
                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [rg.Units];
                    [Y, group_labels_comb] = combineMetadataIndices(rg, object_group, classification_var, pooling_vals);
                    if ~isempty(classification_val) %New division / pool conditions
                        [Y,group_labels_comb] = poolMetadataValues(clf_group_idx, group_labels_comb, classification_val);
                    end
                elseif level == "Recording"
                    input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                    Y = group_idx;
                else
                    error('Unknown level')
                end
            else
                if level == "Unit"
                    Y = arrayfun(@(x) ones(1,length([object_group{x}(1).Units]))*group_idx(x),1:length(group_idx),'un',0);
                    object_group = cellfun(@(x) [x(1).Units],object_group,'un',0);
                    object_group = [object_group{:}];
                else
                    Y = group_idx;
                end
            end

            if ~isempty(feature_names)
                vars = string(input_table.Properties.VariableNames);
                sel_vars = startsWith(vars, feature_names);
                input_table = input_table(:,sel_vars);
            end

            for k = 1:K_fold
                [Y_train, Y_test, train_idx, test_idx] = rg.cv_split(Y, cv, k);

                [X_train, X_test, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);

                [clf,train_acc] = rg.create_classifier(X_train, Y_train, alg, N_hyper);
                [Y_pred,scores] = predict(clf, X_test);

                options = statset('UseParallel',true);
                predImp = oobPermutedPredictorImportance(clf,'Options',options);
% predImp = [];

                result(k).Mdl = clf;
                result(k).Y_pred = Y_pred;
                result(k).Y_test = Y_test;
                result(k).scores = scores;
                result(k).objects = object_group(test_idx);
                result(k).train_acc = train_acc;
                result(k).predImp = predImp;
                result(k).GroupLabels = group_labels_comb;
                result(k).Features = feature_names;
            end
            rg.Classification.(classification_var) = result;
end
