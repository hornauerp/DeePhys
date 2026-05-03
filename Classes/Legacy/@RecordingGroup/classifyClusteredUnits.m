function result = classifyClusteredUnits(rg, true_cluster_idx, feature_groups, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance)
            arguments
                rg RecordingGroup
                true_cluster_idx double
                feature_groups string
                grouping_var string = [] %Metadata field that groups recordings to cultures
                grouping_values = nan %Selected values corresponding to grouping_var
                feature_names string = []
                normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
                normalization_var string = "PlatingDate" % normalization by each value of normalization_var
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
                tolerance double = 1
            end
            network_features = [];
            useClustered = false;
            level = "Unit";
            alg = "rf";

            if isempty(grouping_var)
                object_group = [rg.Recordings];
                input_table = object_group.getUnitFeatures(rg.Units,feature_groups);
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, feature_groups, normalization, useClustered);

            end

            cv = cvpartition(true_cluster_idx,'KFold',K_fold);

            if ~isempty(feature_names)
                vars = string(input_table.Properties.VariableNames);
                sel_vars = startsWith(vars, feature_names);
                input_table = input_table(:,sel_vars);
            end

            for k = 1:K_fold
                [Y_train, Y_test, train_idx, test_idx] = rg.cv_split(true_cluster_idx, cv, k);

                % Balance class distribution by oversampling
                % if oversample
                %     [input_table, Y_train, train_idx, test_idx] = oversample_smaller_class(input_table, Y_train, train_idx);
                % end
                %
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
                result(k).train_acc = train_acc;
                result(k).predImp = predImp;
                result(k).Features = feature_names;
            end
            rg.Classification.UnitClusters = result;
end
