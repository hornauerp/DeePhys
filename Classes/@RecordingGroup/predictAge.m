function result = predictAge(rg, level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, ...
                useClustered, normalization_var, N_hyper, K_fold)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                alg string = "rf" %rf, svm, cnb, knn
                stratification_var = "Mutation" %Specify the variable by which to evenly split training and test dataset
                stratification_values = [] %Corresponding to stratification_var, if a value is specified then this will be used as training data (e.g. train on wildtype, test on mutation)
                pooling_vals cell = {}
                network_features string = "all"
                unit_features string = "all"
                useClustered logical = false
                normalization_var string = "PlatingDate"
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
            end

            t_start = tic;

            if any(arrayfun(@(x) ~isfield(rg.Recordings(x).Metadata,"DIV"),1:length(rg.Recordings)))
                error('Age data missing')
            end

            object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);

            if isempty(stratification_values)
                stratification_var = [stratification_var "DIV"];
                [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var);
                if K_fold == -1 %set kfold to LOOCV
                    K_fold = length(rec_group_idx);
                end
                cv = cvpartition(rec_group_idx,'KFold',K_fold);
                t = templateTree('Surrogate','on','MinLeafSize',5,'NumVariablesToSample','auto');

                for k = 1:K_fold
                    if level == "Unit"
                        train_recs = object_group(cv.training(k));
                        test_recs  = object_group(cv.test(k));
                        train_table = train_recs.getUnitFeatures(unit_features);
                        test_table = test_recs.getUnitFeatures(unit_features);
                        train_idx = logical([ones(1, size(train_table,1)) zeros(1, size(test_table,1))]);
                        test_idx = ~train_idx;
                        input_group = [train_recs.Units, test_recs.Units];
                        input_table = [train_table;test_table];
                        % Build per-unit age vector by expanding each recording's DIV
                        Y_train = arrayfun(@(r) repmat(r.Metadata.DIV, 1, length(r.Units)), train_recs, 'un', 0);
                        Y_train = [Y_train{:}];
                        Y_test  = arrayfun(@(r) repmat(r.Metadata.DIV, 1, length(r.Units)), test_recs, 'un', 0);
                        Y_test  = [Y_test{:}];
                    else
                        input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                        train_idx = cv.training(k);
                        test_idx = cv.test(k);
                        input_group = object_group;
                        metadata = [object_group.Metadata];
                        true_age = [metadata.DIV];
                        Y_train = true_age(train_idx);
                        Y_test = true_age(test_idx);
                    end

                    [X_train, X_test, feature_names] = prepareInputMatrix(rg, input_table, input_group, normalization_var, train_idx, test_idx);

                    clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t);
                    Y_pred = predict(clf, X_test);

                    options = statset('UseParallel',true);
                    predImp = clf.oobPermutedPredictorImportance('Options',options);
                    %                     predImp = [];

                    result(k) = RegressionResult( ...
                        'Mdl', clf, ...
                        'Y_pred', Y_pred, ...
                        'Y_test', Y_test(:), ...
                        'mse_train', resubLoss(clf), ...
                        'objects', object_group(test_idx), ...
                        'predImp', predImp, ...
                        'Features', feature_names, ...
                        'Parameters', struct('algorithm', 'rf', 'regression_var', 'DIV', ...
                            'normalization_var', normalization_var, 'K_fold', K_fold));
                end
            else
                train_group = find(group_labels_comb == stratification_values);

                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [object_group.Units];
                    [unit_group_idx, ~] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);
                    train_idx = unit_group_idx == train_group;
                    test_idx = unit_group_idx ~= train_group;
                    unit_recordings = [object_group.MEArecording];
                    metadata = [unit_recordings.Metadata];
                    true_age = [metadata.DIV];

                elseif level == "Recording"
                    input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                    train_idx = rec_group_idx == train_group;
                    test_idx = rec_group_idx ~= train_group;
                    metadata = [object_group.Metadata];
                    true_age = [metadata.DIV];

                end
                [X_train, X_test] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);
                Y_train = true_age(train_idx);
                Y_test = true_age(test_idx);
                t = templateTree('Surrogate','on','MinLeafSize',5,'NumVariablesToSample','auto');
                clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t,'Options',statset('UseParallel',true));
                Y_pred = predict(clf,X_test);

                options = statset('UseParallel',true);
                predImp = clf.oobPermutedPredictorImportance('Options',options);
%                 predImp = [];

                result = RegressionResult( ...
                    'Mdl', clf, ...
                    'Y_pred', Y_pred, ...
                    'Y_test', Y_test(:), ...
                    'mse_train', resubLoss(clf), ...
                    'objects', object_group(test_idx), ...
                    'predImp', predImp, ...
                    'Features', [], ...
                    'Parameters', struct('algorithm', 'rf', 'regression_var', 'DIV', ...
                        'normalization_var', normalization_var));
            end

            rg.Regression.DIV = result;

            % Log the analysis
            elapsed = toc(t_start);
            summary_metrics = RegressionResult.summarizeFolds(result);
            AnalysisLog.instance().add('predictAge', ...
                struct('level', level, 'K_fold', K_fold), ...
                sprintf('R2=%.3f, Corr=%.3f', summary_metrics.mean_R2, summary_metrics.mean_Correlation), ...
                elapsed);
end
