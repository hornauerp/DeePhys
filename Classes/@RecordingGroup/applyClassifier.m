function result = applyClassifier(rg, level, alg, classification_var, clf_pooling_vals, test_var, test_pooling_vals, network_features, unit_features, feature_names, useClustered,...
                                    grouping_var, grouping_values, normalization, normalization_var, N_hyper, tolerance)
           arguments
               rg RecordingGroup
               level string = "Recording" %Unit or Recording
               alg string = "rf" %rf, svm,
               classification_var string = "Mutation" %Metadata field that determines y_test/y_train
               clf_pooling_vals cell = {}
               test_var string = "Treatment"
               test_pooling_vals cell = {} %First cell corresponds to training data, second to test data
               network_features string = "all"
               unit_features string = "all"
               feature_names string = []
               useClustered logical = false
               grouping_var string = [] %Metadata field that groups recordings to cultures
               grouping_values = nan %Selected values corresponding to grouping_var
               normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
               normalization_var string = "PlatingDate" % normalization by each value of normalization_var
               N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
               tolerance double = 1
           end
           if isempty(grouping_var)
               if level == "Unit"
                   object_group = [rg.Units];
                   input_table = object_group.getUnitFeatures(unit_features);
               elseif level == "Recording"
                   object_group = [rg.Recordings];
                   input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
               else
                   error('Unknown level')
               end
           else
               [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
               if level == "Unit"
                   object_group = cellfun(@(x) [x(1).Units],object_group,'un',0);
                   object_group = [object_group{:}];
               elseif level ~= "Recording"
                   error('Unknown level')
               end
           end

           if ~isempty(feature_names)
               vars = string(input_table.Properties.VariableNames);
               sel_vars = startsWith(vars, feature_names);
               input_table = input_table(:,sel_vars);
           end

           % Split by classification variable
            [Y, group_labels_comb] = combineMetadataIndices(rg, object_group, classification_var, clf_pooling_vals);

            % Identify training and test set
            [test_group_idx, test_group_labels] = combineMetadataIndices(rg, object_group, test_var, test_pooling_vals);
%             [test_group_idx,test_group_labels] = rg.poolMetadataValues(test_group_idx, test_group_labels, );

            train_idx = test_group_idx == 1;
            test_idx = ~train_idx;
            Y_train = Y(train_idx);
            Y_test = Y(test_idx);

            [X_train, X_test] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);
            [clf,train_acc] = rg.create_classifier(X_train, Y_train, alg, N_hyper);
            [Y_pred,scores] = predict(clf, X_test);
            result.Mdl = clf;
            result.Y_pred = Y_pred;
            result.Y_test = Y_test;
            result.scores = scores;
            result.objects = object_group(test_idx);
            result.train_acc = train_acc;
            result.GroupLabels = group_labels_comb;

            rg.Classification.(classification_var).(test_var) = result;
end

function accuracy_mat = assessAppliedClassifier(rg, result, assessment_var, pooling_vals)
            arguments
                rg RecordingGroup
                result struct %Result of applyClassifier function
                assessment_var string %Metadata variable (e.g. "Treatment")
                pooling_vals cell = {}
            end

            [Y_assess, group_labels_comb] = rg.combineMetadataIndices(result.objects, assessment_var, pooling_vals);

            accuracy_mat = nan(length(result.GroupLabels),length(group_labels_comb));
            for t = 1:length(result.GroupLabels)
                for a = 1:length(group_labels_comb)
                    comp_idx = (Y_assess == a) & (result.Y_test == t);
                    accuracy_mat(t,a) = sum(result.Y_test(comp_idx) == result.Y_pred(comp_idx)) ./ sum(comp_idx);
                end
            end
            figure('Color','w');
            heatmap(accuracy_mat,'ColorLimits',[0 1],'XDisplayLabels',group_labels_comb,'YDisplayLabels',result.GroupLabels)
end
