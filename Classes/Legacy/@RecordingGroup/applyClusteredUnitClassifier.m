function result = applyClusteredUnitClassifier(rg, true_cluster_idx, train_idx, test_idx, feature_groups, clf,grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance)
            arguments
               rg RecordingGroup
               true_cluster_idx double %Only needed if clf is not provided
               train_idx logical
               test_idx logical
               feature_groups string
               clf = []
               grouping_var string = [] %Metadata field that groups recordings to cultures
               grouping_values = nan %Selected values corresponding to grouping_var
               feature_names string = []
               normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
               normalization_var string = "PlatingDate" % normalization by each value of normalization_var
               N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
               tolerance double = 0
            end
            network_features = [];
            useClustered = false;
            level = "Unit";
            alg = "rf";

            if isempty(grouping_var)
                object_group = [rg.Recordings];
                input_table = object_group.getUnitFeatures(rg.Units, feature_groups);
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, feature_groups, normalization, useClustered);

            end

            if ~isempty(feature_names)
                vars = string(input_table.Properties.VariableNames);
                sel_vars = matches(vars, feature_names);
                input_table = input_table(:,sel_vars);
            end

            [X_train, X_test, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);

            if isempty(clf)
                [clf,train_acc] = rg.create_classifier(X_train, true_cluster_idx(train_idx), alg, N_hyper);
                result.train_acc = train_acc;
            end

            [Y_pred,scores] = predict(clf, X_test);

            result.Y_pred = Y_pred;
            result.scores = scores;
            result.feature_names = feature_names;
            result.Mdl = clf;
end
% assignUnitClusterIdx: moved to assignUnitClusterIdx.m
