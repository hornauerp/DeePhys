function result = classifyByFeatureGroups(rg, level, alg, classification_var, pooling_vals, network_features, unit_features, feature_names, useClustered, ...
                                    grouping_var, grouping_values, normalization, normalization_var, N_hyper, K_fold, tolerance, cv_grouping)
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
                cv_grouping string = "recording" % "none" (standard), "recording", or "culture" — GroupedCV split level
            end

            t_start = tic;

            if isempty(grouping_var) %Get group indices for stratification
                object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            end
            [group_idx, group_labels_comb] = rg.combineMetadataIndices(object_group, classification_var, pooling_vals);


            if isempty(grouping_var)
                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [rg.Units];
                    [Y, group_labels_comb] = combineMetadataIndices(rg, object_group, classification_var, pooling_vals);
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

            % Build CV partition — hierarchy-aware or standard
            use_grouped = cv_grouping ~= "none";
            if use_grouped
                if iscell(Y)
                    Y_flat = [Y{:}];
                else
                    Y_flat = Y;
                end
                switch cv_grouping
                    case "recording"
                        gcv = GroupedCV.byRecording(object_group, K_fold, Y_flat);
                    case "culture"
                        gcv = GroupedCV.byCulture(object_group, K_fold, Y_flat);
                    otherwise
                        error('Unknown cv_grouping: "%s". Use "none", "recording", or "culture".', cv_grouping);
                end
                K_fold = gcv.NumFolds;
            else
                if K_fold == -1
                    K_fold = length(group_idx);
                end
                cv = cvpartition(group_idx, 'KFold', K_fold);
            end

            for k = 1:K_fold
                if use_grouped
                    [train_idx, test_idx] = gcv.fold(k);
                    if iscell(Y)
                        Y_train = [Y{train_idx}]; Y_test = [Y{test_idx}];
                    else
                        Y_train = Y(train_idx); Y_test = Y(test_idx);
                    end
                else
                    [Y_train, Y_test, train_idx, test_idx] = rg.cv_split(Y, cv, k);
                end

                [X_train, X_test, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);

                [clf,train_acc] = rg.create_classifier(X_train, Y_train, alg, N_hyper);
                [Y_pred,scores] = predict(clf, X_test);

                if alg == "rf"
                    options = statset('UseParallel',true);
                    predImp = oobPermutedPredictorImportance(clf,'Options',options);
                else
                    predImp = [];
                end
% predImp = [];

                result(k) = ClassificationResult( ...
                    'Mdl', clf, ...
                    'Y_pred', Y_pred, ...
                    'Y_test', Y_test, ...
                    'scores', scores, ...
                    'objects', object_group(test_idx), ...
                    'train_acc', train_acc, ...
                    'predImp', predImp, ...
                    'GroupLabels', group_labels_comb, ...
                    'Features', feature_names, ...
                    'Parameters', struct('algorithm', alg, 'classification_var', classification_var, ...
                        'normalization_var', normalization_var, 'K_fold', K_fold, 'N_hyper', N_hyper, ...
                        'cv_grouping', cv_grouping));
            end
            rg.Classification.(classification_var) = result;

            % Log the analysis
            elapsed = toc(t_start);
            summary_metrics = ClassificationResult.summarizeFolds(result);
            AnalysisLog.instance().add('classifyByFeatureGroups', ...
                struct('level', level, 'alg', alg, 'classification_var', classification_var, ...
                    'K_fold', K_fold, 'cv_grouping', cv_grouping), ...
                sprintf('Acc=%.3f+/-%.3f, F1=%.3f', summary_metrics.mean_accuracy, ...
                    summary_metrics.std_accuracy, summary_metrics.mean_F1), ...
                elapsed);
end
