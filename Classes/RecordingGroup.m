classdef RecordingGroup < handle
    
    %%%%%%%TODO%%%%%%%
    % 
    %%%%%%%%%%%%%%%%%%
    
    properties
        Recordings
        Units
        Cultures
        Parameters %Inclusion // exclusion criteria
        DimensionalityReduction
        Clustering
        Classification        
    end
    
    methods (Hidden)
        function rg = RecordingGroup(recording_array, parameters)
            arguments
                recording_array MEArecording
                parameters struct = struct();
            end
            
            rg.parseParameters(parameters);
            keep_recordings = rg.filterRecordingArray([recording_array.Metadata]);
            rg.Recordings = recording_array(keep_recordings);
            rg.Units = [rg.Recordings.Units];
            if isfield([rg.Recordings.Metadata],'ChipID')
                rg.Cultures = rg.groupCultures();
            end
            fprintf('Initialized RecordingGroup with %i recordings and %i units\n',...
                length(rg.Recordings),length(rg.Units))
        end
        
        function parseParameters(rg, parameters)
            rg.Parameters = rg.returnDefaultParams();
            parameter_fields = fieldnames(parameters);
            if isempty(parameter_fields)
                warning("No parameters provided, using default ones")
            else
                for pf = 1:length(parameter_fields)
                    if ismember(parameter_fields{pf},fieldnames(rg.Parameters))
                        subfields = fieldnames(parameters.(parameter_fields{pf}));
                        for sf = 1:length(subfields)
                            if ismember(subfields{sf},fieldnames(rg.Parameters.(parameter_fields{pf})))
                                rg.Parameters.(parameter_fields{pf}).(subfields{sf}) = parameters.(parameter_fields{pf}).(subfields{sf});
                            else
                                warning("Unrecognized parameter field: %s.%s",parameter_fields{pf}{:},subfields{sf}{:})
                            end
                        end
                    else
                        warning("Unrecognized parameter field: %s",pf{:})
                    end
                end
                disp("Imported custom parameters")
            end
        end
        
        function filtered_idx = filterRecordingArray(rg, metadata_array, inclusion, exclusion)
            arguments
                rg RecordingGroup
                metadata_array
                inclusion = rg.Parameters.Selection.Inclusion;
                exclusion = rg.Parameters.Selection.Exclusion;
            end
            
            % Input as cell arrays
            for i = 1:length(inclusion)
                if isnumeric(inclusion{i}{2}) || isstring((inclusion{i}{2}))
                    value_vector = [metadata_array.(inclusion{i}{1})];
                    idx = ismember(value_vector,[inclusion{i}{2:end}]);
                elseif iscell(inclusion{i}{2})
                    value_vector = {metadata_array.(inclusion{i}{1})};
                    idx = ismember(value_vector,inclusion{i}{2:end});
                else
                    value_vector = [metadata_array.(inclusion{i}{1})];
                    idx = ismember(value_vector,inclusion{i}(2:end));
                end
                if ~exist('inclusion_idx','var')
                    inclusion_idx = idx;
                else
                    inclusion_idx = idx & inclusion_idx;
                end
            end
            if ~exist('inclusion_idx','var')
               inclusion_idx = ones(1,length(metadata_array)); 
            end
            
            for i = 1:length(exclusion)
                if isnumeric(exclusion{i}{2}) || isstring((exclusion{i}{2}))
                    value_vector = [metadata_array.(exclusion{i}{1})];
                    idx = ismember(value_vector,[exclusion{i}{2:end}]);
                else
                    value_vector = {metadata_array.(exclusion{i}{1})};
                    idx = ismember(value_vector,exclusion{i}(2:end));
                end
                if ~exist('exclusion_idx','var')
                    exclusion_idx = idx;
                else
                    exclusion_idx = idx | exclusion_idx;
                end
            end
            if ~exist('exclusion_idx','var')
               exclusion_idx = zeros(1,length(metadata_array)); 
            end
            
            filtered_idx = inclusion_idx & ~exclusion_idx;
        end
        
        function cultures = groupCultures(rg)
            cultures = [];
            metadata = [rg.Recordings.Metadata];
            if ~isfield(metadata,"ChipID")
                warning('No chip IDs provided, cannot group recordings')
                return
            end
            chip_ids = unique([metadata.ChipID]);
            if isempty(chip_ids)
                warning('No chip IDs provided, cannot group recordings')
                return
            end
            for id = chip_ids
                inclusion = {{'ChipID',id}};
                chip_idx = rg.filterRecordingArray(metadata, inclusion);
                chip_array = rg.Recordings(chip_idx);
                chip_metadata = [chip_array.Metadata];
                if ~isfield(metadata,"PlatingDate")
                    warning('No plating date provided, cannot differentiate batches')
                    return
                end
                chip_plating_dates = unique([chip_metadata.PlatingDate]);
                for pd = chip_plating_dates
                   inclusion = {{'PlatingDate',pd}};
                   plating_idx = rg.filterRecordingArray(chip_metadata, inclusion);
                   culture_array = chip_array(plating_idx);
                   cultures = [cultures {culture_array}];
                end
            end
        end
    end
    
    methods (Static)
       function defaultParams = returnDefaultParams()
           %Input to the filterRecordingArray function
           defaultParams.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
           defaultParams.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
       end
       
       function aligned_wf = alignWaveforms(unit_array)
           ref_wf = [unit_array.ReferenceWaveform];
           ref_wf = ref_wf(sum(ref_wf,2)~=0,:);
           [~,i] = min(ref_wf,[],1);
           peak_idx = mean(i);
           max_offset = round(peak_idx/2);
           x = max_offset:size(ref_wf,1)+max_offset-1;
           xq = 1:size(ref_wf,1)+2*max_offset;
           
           interp_wf = interp1(x,ref_wf,xq,"linear",'extrap');
           rm_idx = find(ceil(abs(peak_idx - i)) >= max_offset);
           aligned_wf = zeros(size(ref_wf));
           for r = 1:size(interp_wf,2)
               if ~any(rm_idx == r)
                   start_idx = i(r) - max_offset;
                   aligned_wf(:,r) = interp_wf(start_idx:start_idx + size(ref_wf,1) - 1 ,r);
               else
                   aligned_wf(:,r) = ref_wf(:,r);
               end
           end
       end
       
       function [clf,train_acc] = create_classifier(X_train,Y_train,alg,N_hyper)
           if N_hyper > 0
               switch alg
                   case 'svm'
                       clf = fitcsvm(X_train,Y_train,'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
                           struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                           'Verbose',0));
                   case 'cnb'
                       clf = fitcnb(X_train,Y_train,'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
                           struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                           'Verbose',0));
                   case 'knn'
                       clf = fitcknn(X_train,Y_train,'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
                           struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                           'Verbose',0));
                   case 'rf'
                       hyperparams = {'NumLearningCycles','MinLeafSize','MaxNumSplits','SplitCriterion','NumVariablesToSample'};
                       t = templateTree('Reproducible',true);
                       clf = fitcensemble(X_train,Y_train,'Method','Bag','OptimizeHyperparameters',hyperparams,'Learners',t, ...
                           'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                           'Verbose',0));
               end
               train_acc = clf.HyperparameterOptimizationResults.MinObjective;
           else
               switch alg
                   case 'svm'
                       clf = fitcsvm(X_train,Y_train);
                   case 'cnb'
                       clf = fitcnb(X_train,Y_train);
                   case 'knn'
                       clf = fitcknn(X_train,Y_train);
                   case 'rf'
                       t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all','Reproducible',true);
                       clf = fitcensemble(X_train,Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t,'Options',statset("UseParallel",true));
               end
               train_acc = 1-resubLoss(clf,'LossFun','classiferror');
           end
       end
       
       function [Y_train, Y_test, train_idx, test_idx] = cv_split(Y, cv, k)
          arguments
             Y
             cv cvpartition
             k double
          end
          
           if iscell(Y)
               Y_train = [Y{cv.training(k)}];
               Y_test = [Y{cv.test(k)}];
               k_train = cv.training(k);
               train_idx = arrayfun(@(x) ones(size(Y{x})) * k_train(x),1:length(Y),'un',0);
               train_idx = [train_idx{:}];
               test_idx = ~train_idx;
               
           else
               Y_train = Y(cv.training(k));
               Y_test = Y(cv.test(k));
               train_idx = cv.training(k);
               test_idx = cv.test(k);
           end
       end
       
       function [new_group_idx,new_group_labels] = poolMetadataValues(group_idx, group_labels, classification_val)
           arguments
               group_idx double %Original group idx, corresponding to group_labels
               group_labels string
               classification_val string %Value(s) that are to be pooled (against) -> classification_val vs rest
           end
           clf_group_idx = find(ismember(group_labels, classification_val));
           new_group_idx = (group_idx == clf_group_idx) * 1;
           new_group_labels(1) = join(group_labels(clf_group_idx));
           clf_group_idx = find(~ismember(group_labels, classification_val));
           new_group_idx(new_group_idx == 0) = 2;
           new_group_labels(2) = join(group_labels(clf_group_idx));
       end
    end
    
    methods 
        
        function [feature_table, culture_array] = aggregateCultureFeatureTables(rg, level, metadata_field, values, tolerance, network_features, unit_features, normalization, useClustered)
            arguments
                rg RecordingGroup
                level string = "Recording"
                metadata_field string = "DIV"
                values {isnumeric} = [7,14,21,28] %refers to values of metadata_field // NaN to use the maximum number of available timepoints
                tolerance {isnumeric} = 1 %deviation from the actual age that will still be considered valid (e.g. age = 7 and tolerance = 1 considers DIVs 6-8)
                network_features string = "all"
                unit_features string = "all"
                normalization string = "baseline" % "baseline" (divided by first timepoint) or "scaled (scaled between 0 and 1)
                useClustered logical = false
            end
            
            rec_metadata = [rg.Recordings.Metadata];
            unique_values = unique([rec_metadata.(metadata_field)]);
            if isnan(values)
                values = unique_values;
            end
            value_dev = abs(values - unique_values');
            value_count = sum(value_dev<=tolerance);
            if any(value_count == 0)
               warning("DIV " + num2str(values(value_count == 0)) + " without corresponding cultures")
            end
            values = values(value_count > 0);
            if isempty(values)
               error('No recordings found') 
            end
            culture_table_array = {};
            culture_array = {};
            for iC = 1:length(rg.Cultures)
                rec_table_array = cell(1,length(values));
                culture_metadata = [rg.Cultures{iC}.Metadata];
                value = [culture_metadata.(metadata_field)];
                if length(value) >= length(values)
                    value_dev = abs(value - values');
                    value_count = sum(value_dev<=tolerance,1);
                    value = value(value_count==1);
                    if sum(value_count==1) == length(values) % Need to find a way to handle two recordings falling within the tolerance window
                        sel_rec = rg.Cultures{iC}(value_count == 1);
                        culture_array = [culture_array {sel_rec}];
                        for iR = 1:length(sel_rec)
                            if level == "Unit"
                                iR_table = getUnitFeatures(rg.Cultures{iC}(iR), unit_features);
                            elseif level == "Recording"
                                iR_table = getRecordingFeatures(rg.Cultures{iC}(iR), network_features, unit_features, useClustered);
                            else 
                                error('Unknown level, select either "Unit" or "Recording"')
                            end
                            iR_table.Properties.VariableNames = iR_table.Properties.VariableNames + "_" + string(value(iR));
                            rec_table_array{iR} = iR_table;
                        end
                        if normalization == "baseline"
                            norm_mat = arrayfun(@(x) rec_table_array{x}.Variables./rec_table_array{1}.Variables,1:length(rec_table_array),'un',0);
                            norm_mat = [norm_mat{2:end}];
                            norm_mat(isnan(norm_mat)) = 0;
                            norm_mat(isinf(norm_mat)) = max(norm_mat(norm_mat < Inf),[],'all');
                            
                            norm_table = [rec_table_array{2:end}];
                            norm_vars = norm_table.Properties.VariableNames;
                            culture_table_array{iC} = array2table(norm_mat,'VariableNames',norm_vars);
                            
                        elseif normalization == "scaled"
                            norm_mat = cellfun(@(x) x.Variables,rec_table_array,'un',0);
                            norm_mat = cat(3,norm_mat{:});
                            norm_mat = normalize(norm_mat,3,'range',[0 1]);
                            re_mat = reshape(norm_mat,size(norm_mat,1),[]);
                            norm_table = [rec_table_array{:}];
                            norm_vars = norm_table.Properties.VariableNames;
                            culture_table_array{iC} = array2table(re_mat,'VariableNames',norm_vars);
                        else
                            culture_table_array{iC} = [rec_table_array{:}];
                        end
                    else
                        continue
                    end
                else
                   continue
               end
            end
%             clean_culture_tables = ~cellfun(@isempty, culture_table_array);
            feature_table = vertcat(culture_table_array{:});
            keep_idx = ~isnan(std(feature_table.Variables,'omitnan')); %Remove variables with 0 variance
            feature_table = feature_table(:,keep_idx);
        end
        
        function [norm_train_data, norm_test_data] = normalizeByGroup(rg, feat_mat, object_group, grouping_var, train_idx, test_idx)
            arguments
                rg RecordingGroup
                feat_mat {isnumeric}
                object_group %Units/recordings/cultures corresponding to the rows of feat_mat
                grouping_var string = "PlatingDate" %Has to correspond to a MEArecording metadata field
                train_idx (1,:) logical = logical(ones(1,size(feat_mat,1))) %Default to unsupervised // should be binary and length of object_group
                test_idx (1,:) logical = logical(zeros(1,size(feat_mat,1)))
            end
            
            train_idx = logical(train_idx);
            test_idx = logical(test_idx);
            
            switch class(object_group)
                case 'Unit'
                    recordings = [object_group.MEArecording];
                    metadata = [recordings.Metadata];
                    metadata = [metadata.(grouping_var)];
                    
                case 'MEArecording'
                    metadata = [object_group.Metadata];
                    metadata = [metadata.(grouping_var)];
                    
                case 'cell'
                    metadata = cellfun(@(x) string(x(1).Metadata.(grouping_var)),object_group);
            end
            
            
            [iG,G] = findgroups(metadata);
            for g = unique(iG)
                iBatch_train = (iG == g);% & train_idx;
                [feat_mat(iBatch_train,:), batch_mean, batch_sd] = normalize(feat_mat(iBatch_train,:));
%                 if any(test_idx)
%                     iBatch_test = iG == g & test_idx;
%                     feat_mat(iBatch_test,:) = normalize(feat_mat(iBatch_test,:),'center',batch_mean,'scale',batch_sd);
%                 end
            end
            if length(G) > 1 %Only normalize again if more than 1 group exists
                [norm_train_data, train_mean, train_sd] = normalize(feat_mat(train_idx,:));
                norm_test_data = normalize(feat_mat(test_idx,:),'center',train_mean,'scale',train_sd);
            else
                norm_train_data = feat_mat(train_idx,:);
                norm_test_data = feat_mat(test_idx,:);
            end
        end
        
        function [norm_train, norm_test] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx)
            arguments
                rg RecordingGroup
                input_table table
                object_group %Units/recordings/cultures corresponding to the rows of feat_mat
                normalization_var string = "PlatingDate" %Has to correspond to a MEArecording metadata field
                train_idx (1,:) logical = logical(ones(1,size(input_table,1))) %Default to unsupervised // should be logical and length of object_group
                test_idx (1,:) logical = logical(zeros(1,size(input_table,1)))
            end
            
            input_mat = input_table.Variables;
            input_mat(isnan(input_mat)) = 0; %Handle rare case where NaN appears
            if isempty(normalization_var)
                [norm_train, mean_train, sd_train] = normalize(input_mat(train_idx,:));
                nan_idx = any(isnan(norm_train));
                norm_train(:,nan_idx) = [];
                norm_train = norm_train./max(abs(norm_train)); %Scale data
                if sum(test_idx) > 0
                    norm_test = normalize(input_mat(test_idx,:),'center',mean_train,'scale',sd_train);
                    norm_test(:,nan_idx) = [];
                    norm_test = norm_test./max(abs(norm_train));
                else
                    norm_test = [];
                end
            else
                [norm_train, norm_test] = normalizeByGroup(rg, input_mat, object_group, normalization_var, train_idx, test_idx); %Normalize
                nan_idx = any(isnan(norm_train)) | any(isnan(norm_test),1);
                norm_train(:,nan_idx) = [];
                
                norm_train = norm_train./max(abs(norm_train)); %Scale data
                if sum(test_idx) > 0
                    norm_test(:,nan_idx) = [];%Remove peak NaNs
                    norm_test = norm_test./max(abs(norm_train)); %Scale data
                else
                    norm_test = [];
                end
            end
            
        end
        
        function reduction = reduceDimensionality(rg, level, method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance)
            arguments
               rg RecordingGroup
               level string = "Unit" %Unit, Recording or Culture level
               method string = "UMAP" %UMAP, PCA
               n_dims (1,1) {isnumeric} = 2 %Number of output dimensions
               unit_features string = ["ReferenceWaveform","ActivityFeatures"]
               network_features string = "all"
               useClustered logical = false
               normalization_var string = "PlatingDate"
               grouping_var string = [] %Has to correspond to a MEArecording metadata field
               grouping_values {isnumeric} = [7,14,21,28] %Only relevant for DimensionalityReduction on the culture level
               normalization string = "scaled" %"scaled" or "baseline" or []
               tolerance {isnumeric} = 1 %Gives tolerance for culture selection by age (e.g. age=7 tolerance=1 allows DIVs 6-8)
            end
            
%             fprintf('Performing %s dimensionality reduction on %ss and normalizing by %s\n', method, level, grouping_var)
            if isempty(grouping_var) %Get group indices for stratification
                object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            end
            
            if isempty(grouping_var)
                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [rg.Units];
                elseif level == "Recording"
                    input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                else 
                    error('Unknown level')
                end
            else
                if level == "Unit"
                    object_group = cellfun(@(x) [x(1).Units],object_group,'un',0);
                    object_group = [object_group{:}];
                end
            end
            
            norm_data = prepareInputMatrix(rg, input_table, object_group, normalization_var);
            
            switch method
                case "UMAP"
                    color_file = fullfile(rg.Recordings(1).getParentPath(),'umap','colorsByName.properties');
                    [reduction, umap, clusterIdentifiers, extras] = run_umap(norm_data,'n_components',n_dims,'n_neighbors',100,'min_dist',0.1,'cluster_detail','adaptive','spread',1,'sgd_tasks',20,...
                        'verbose','none','color_file',color_file);
%                     [reduction, umap, clusterIdentifiers, extras] = run_umap(norm_data,'n_components',n_dims,'sgd_tasks',20,...
%                         'verbose','none','color_file',color_file);
                case "PCA"
                    [coeff,reduction,latent] = pca(norm_data);
                    
                case "tSNE"
                    reduction = tsne(norm_data,'Algorithm','exact','Exaggeration',10,'Perplexity',30);
            end
            reduction = reduction(:,1:n_dims);
            
            rg.DimensionalityReduction.(level).(method).Reduction = reduction;
            rg.DimensionalityReduction.(level).(method).GroupingVariable = grouping_var;
            rg.DimensionalityReduction.(level).(method).UnitFeatures = unit_features;
            rg.DimensionalityReduction.(level).(method).NetworkFeatures = network_features;
            rg.DimensionalityReduction.(level).ObjectGroup = object_group;
%             rg.plot_dimensionality_reduction(reduction);
        end
        
        function result = predictAge(rg, level, alg, stratification_var, stratification_values, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                alg string = "rf" %rf, svm, cnb, knn
                stratification_var = "Mutation" %Specify the variable by which to split training and test dataset (e.g. train on wildtype, test on mutation)
                stratification_values = [] %Corresponding to stratification_var, if a value is specified then this will be used as training data
                network_features string = "all"
                unit_features string = "all"
                useClustered logical = false
                normalization_var string = "PlatingDate"
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
            end
            
            if any(arrayfun(@(x) ~isfield(rg.Recordings(x).Metadata,"DIV"),1:length(rg.Recordings)))
                error('Age data missing')
            end
            
            object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var);
            
            if isempty(stratification_values)
                stratification_var = [stratification_var "DIV"];
                [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var);
                cv = cvpartition(rec_group_idx,'KFold',K_fold);
                t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');
                
                for k = 1:K_fold
                    if level == "Unit"
                        train_table = object_group(cv.training(k)).getUnitFeatures(unit_features);
                        test_table = object_group(cv.test(k)).getUnitFeatures(unit_features);
                        train_idx = logical([ones(1, size(train_table,1)) zeros(1, size(test_table,1))]);
                        test_idx = ~train_idx;
                        input_group = [object_group(cv.training(k)).Units, object_group(cv.test(k)).Units];
                        input_table = [train_table;test_table];
                        unit_recordings = [object_group.MEArecording];
                        metadata = [unit_recordings.Metadata];
                        true_age = [metadata.DIV];
                        Y_train = true_age(train_idx);
                        Y_test = true_age(test_idx);
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
                    
                    [X_train, X_test] = prepareInputMatrix(rg, input_table, input_group, normalization_var, train_idx, test_idx);
                    
                    clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t);
                    Y_pred = predict(clf, X_test);
                    predImp = clf.oobPermutedPredictorImportance;
                    
                    result(k).Mdl = clf;
                    result(k).Y_pred = Y_pred;
                    result(k).Y_test = Y_test;
                    result(k).mse_train = resubLoss(clf);
                    result(k).predImp = predImp;
                    result(k).GroupLabels = group_labels_comb;
                end
            else
                train_group = find(group_labels_comb == stratification_values);
                
                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [object_group.Units];
                    [unit_group_idx, ~] = combineMetadataIndices(rg, object_group, stratification_var);
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
                t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');
                clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t);
                Y_pred = predict(clf,X_test);
                predImp = clf.oobPermutedPredictorImportance;
                
            end
            
            
        end
        
        function result = classifyByFeatureGroups(rg, level, alg, classification_var, classification_val, network_features, unit_features, useClustered,... 
                                    grouping_var, grouping_values, normalization, normalization_var, N_hyper, K_fold, tolerance)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                alg string = "rf" %rf, svm, 
                classification_var string = "Mutation" %Metadata field that determines y_test/y_train
                classification_val string = [] %Values corresponding to classification_var that is being classified for(e.g. "LNA")
                network_features string = "all"
                unit_features string = "all"
                useClustered logical = false
                grouping_var string = [] %Metadata field that groups recordings to cultures
                grouping_values = nan %Selected values corresponding to grouping_var
                normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
                normalization_var string = "PlatingDate" % normalization by each value of normalization_var
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
                tolerance double = 1
            end
            
            if isempty(grouping_var) %Get group indices for stratification
                object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            end
            [group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, classification_var);
            
            if ~isempty(classification_val) %New division / pool conditions
                [group_idx,group_labels_comb] = poolMetadataValues(group_idx, group_labels_comb, classification_val);
            end
            
            if K_fold == -1 %set kfold to LOOCV
                K_fold = length(group_idx);
            end
            cv = cvpartition(group_idx,'KFold',K_fold);
            
            if isempty(grouping_var)
                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [rg.Units];
                    [Y, group_labels_comb] = combineMetadataIndices(rg, object_group, classification_var);
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
                      
            
            for k = 1:K_fold
                [Y_train, Y_test, train_idx, test_idx] = rg.cv_split(Y, cv, k);
                
                [X_train, X_test] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);

%                 feat_mat = input_table.Variables;
%                 feat_mat(isnan(feat_mat)) = 0;
%                 feat_mat = feat_mat ./ max(abs(feat_mat));
%                 norm_mat = normalize(feat_mat);
%                 X_train = norm_mat(train_idx,:);
%                 X_test = norm_mat(test_idx,:);
                
                [clf,train_acc] = rg.create_classifier(X_train, Y_train, alg, N_hyper);
                [Y_pred,scores] = predict(clf, X_test);
                %                 predImp = clf.oobPermutedPredictorImportance();
                predImp = [];
                
                result(k).Mdl = clf;
                result(k).Y_pred = Y_pred;
                result(k).Y_test = Y_test;
                result(k).scores = scores;
                result(k).objects = object_group(test_idx);
                result(k).train_acc = train_acc;
                result(k).predImp = predImp;
                result(k).GroupLabels = group_labels_comb;
            end
            
        end
        
        function result = applyClassifier(rg, level, alg, classification_var, classification_val, test_var, test_val, network_features, unit_features, useClustered,... 
                                    grouping_var, grouping_values, normalization, normalization_var, N_hyper, tolerance)
           arguments
               rg RecordingGroup
               level string = "Recording" %Unit or Recording
               alg string = "rf" %rf, svm,
               classification_var string = "Mutation" %Metadata field that determines y_test/y_train
               classification_val string = [] %Values corresponding to classification_var that are pooled (against) (e.g. "LNA" would classify all others vs LNA)
               test_var string = "Treatment"
               test_val string = "LNA"
               network_features string = "all"
               unit_features string = "all"
               useClustered logical = false
               grouping_var string = [] %Metadata field that groups recordings to cultures
               grouping_values = nan %Selected values corresponding to grouping_var
               normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
               normalization_var string = "PlatingDate" % normalization by each value of normalization_var
               N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
               tolerance double = 1
           end
           
           if isempty(grouping_var) %Get group indices for stratification
                object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            end
            [group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, classification_var);
            if ~isempty(classification_val)
                [group_idx,group_labels_comb] = rg.poolMetadataValues(group_idx, group_labels_comb, classification_val);
            end
            
            [test_group_idx, test_group_labels] = combineMetadataIndices(rg, object_group, test_var);
            if ~isempty(test_val)
                [test_group_idx,test_group_labels] = rg.poolMetadataValues(test_group_idx, test_group_labels, test_val);
            end
        end
        
        function regressionByFeatures(rg)
           arguments
               rg RecordingGroup
           end
        end
        
        function [iMetadata, metadata_groups] = returnMetadataArray(rg, metadata_object, metadata_name)
           arguments
               rg RecordingGroup
               metadata_object %Can be a string referring to the respective group (Unit,Recording,Culture) within the RecordingGroup, or a Unit,Recording or Culture array
               metadata_name string
           end
           
           if isstring(metadata_object)
              switch metadata_object
                  case "Unit"
                      metadata_object = [rg.Units.MEArecording];
                      
                  case "Recording"
                      metadata_object = [rg.Recordings];
                      
                  case "Culture"
                      metadata_object = cellfun(@(x) x(1),rg.Cultures);
                      
              end
              
           elseif iscell(metadata_object)
               metadata_object = cellfun(@(x) x(1),metadata_object);
               
           elseif class(metadata_object) == "Unit"
               metadata_object = [metadata_object.MEArecording];
           end
           
           metadata_struct = [metadata_object.Metadata];
           assert(isfield(metadata_struct,metadata_name),"Not a valid metadata field")
           metadata = [metadata_struct.(metadata_name)];
           [iMetadata, metadata_groups] = findgroups(metadata);
        end
        
        function [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, metadata_names)
            arguments
                rg RecordingGroup
                metadata_object %Can be a string referring to the respective group (Unit,Recording,Culture) within the RecordingGroup, or a Unit,Recording or Culture array
                metadata_names string
            end
            
            idx_array = cell(1,length(metadata_names));
            group_array = cell(1,length(metadata_names));
            for m = 1:length(metadata_names)
                [idx_array{end-m+1}, group_array{end-m+1}] = returnMetadataArray(rg, metadata_object, metadata_names(m));
            end
            group_vals = cellfun(@unique, idx_array, "un",0);
            group_combs = combvec(group_vals{:});
            group_idx = vertcat(idx_array{:});
            [~,cluster_idx] = ismember(group_idx', group_combs','rows');
            group_labels = arrayfun(@(x) group_array{x}(group_combs(x,:))',1:length(group_array),"un",0);
            if length(metadata_names) == 1
                group_labels_comb = [group_labels{:}];
            else
                group_labels_comb = join([group_labels{end:-1:1}]);
                cluster_ids = unique(cluster_idx);
                max_clust = length(group_labels_comb);
                group_labels_comb = group_labels_comb(cluster_ids);
                new_ids = cumsum(sum(cluster_ids == 1:max_clust)); %Remove ids without corresponding group
                cluster_idx = new_ids(cluster_idx);
            end
            if ~isstring(group_labels_comb)
                group_labels_comb = string(group_labels_comb);
            end
        end
        
        function cluster_idx = clusterByFeatures(rg, input_type, level, method, N_clust, grouping_var)
           arguments
                rg RecordingGroup
                input_type string = "UMAP" %"RAW","UMAP","PCA"
                level string = "Unit" %Unit, Recording, Culture
                method string = "kmeans" %kmeans, hierarchical, spectral, gmm
                N_clust {isnumeric} = 0 %Number of imposed clusters // 0 selects N_clust by finding the optimal silhouette score
                grouping_var string = "PlatingDate" % only for input_type "RAW"
           end
           
           if strcmp(input_type, "RAW")
               object_group = [rg.Units];
               unit_feature_array = [object_group.AverageFeatures];
               feature_table = rg.Recordings.concatenateFeatures(unit_feature_array, "all");
               input_mat = feature_table.Variables;
               input_mat(isnan(input_mat)) = 0;%Handle rare case where NaN appears
               norm_data = normalizeByGroup(rg, input_mat, object_group, grouping_var); %Normalize
               feat_mat = norm_data./max(abs(norm_data)); %Scale data
           elseif contains(input_type, ["PCA","UMAP"])
               feat_mat = rg.DimensionalityReduction.(level).(input_type).Reduction;
           else
              error("Unknown input type") 
           end
           
           switch method
               case "kmeans"
                   cluster_fun = @(x,k) kmeans(x,k); 
                   
               case "hierarchical"
                   cluster_fun = @(x,k) clusterdata(x,k); 
                   
               case "dbscan"
                   error("Not yet implemented")
                   
               case "spectral"
                   cluster_fun = @(x,k) spectralcluster(x,k);
                   
               case "gmm"
                   assert(size(feat_mat,1)>size(feat_mat,2),"GMM needs more samples than variables")
                   cluster_fun = @(x,k) cluster(fitgmdist(x,k),x);
                   
           end
           if N_clust == 0
               evaluation = evalclusters(feat_mat,cluster_fun,'silhouette',"KList",1:6);
               N_clust = evaluation.OptimalK;
               cluster_idx = cluster_fun(feat_mat,N_clust);
           else
               cluster_idx = cluster_fun(feat_mat,N_clust);
           end
           rg.Clustering.(level).(method).Index = cluster_idx;
           rg.Clustering.(level).(method).k = N_clust;
           rg.Clustering.(level).(method).Input = input_type;
           
           figure('Color','w');
           plot_cluster_outlines(rg,feat_mat, cluster_idx)
           title(sprintf('%s clustering on %ss',method, level))
        end
        
        
        
        function assignUnitClusterIdx(rg,method,calc_feat)
           arguments
               rg RecordingGroup
               method string = "kmeans"
               calc_feat logical = true %(Re)calculate unit feature averages per cluster
           end
           cluster_idx = rg.Clustering.Unit.(method).Index;
           cluster_idx = num2cell(cluster_idx); %Prepare to use deal to assign cluster ids
           [rg.DimensionalityReduction.Units.ClusterID] = deal(cluster_idx{:});
           
           if calc_feat
               N_clust = num2cell(ones(size(rg.Recordings))*max([cluster_idx{:}]));
               [rg.Recordings.NumUnitClusters] = deal(N_clust{:});
               for r = 1:length(rg.Recordings)
                  rg.Recordings(r).calculateClusterSingleCellFeatures();
               end
           end
        end
        
        function removeUnitsByCluster(rg, method, ID, concatenated)
            arguments
               rg RecordingGroup
               method string
               ID double
               concatenated = false % If true, units are removed for each recording of a culture and not only for individual recordings
            end
            rg.assignUnitClusterIdx(method,false);
            
            if concatenated
                for c = 1:length([rg.Cultures])
                    rm_idx = [];
                    for r = 1:length([rg.Cultures{c}])
                        rm_idx = [rm_idx find([rg.Cultures{c}(r).Units.ClusterID] == ID)];
                    end
                    rm_idx = unique(rm_idx);
                    for r = 1:length([rg.Cultures{c}])
                        rg.Cultures{c}(r).Units(rm_idx) = [];
                    end
                end
            else
                
                for r = 1:length([rg.Recordings])
                    units = [rg.Recordings(r).Units];
                    rg.Recordings(r).Units([units.ClusterID] == ID) = [];
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [cluster_idx, group_labels_comb] = plot_true_clusters(rg, level, method, grouping_var, nodeSz, mapSz, sigma)
            arguments
                rg RecordingGroup
                level string = "Unit" %Unit, Recording, Culture
                method string = "UMAP" %dim reduction method: "UMAP","PCA"
                grouping_var string = "Mutation"
                nodeSz (1,1) {isnumeric} = 10
                mapSz {isnumeric} = 300
                sigma (1,1) {isnumeric} = mapSz/60
%                 cmap (:,3) {isnumeric} = othercolor('Set19',max(cluster_idx));
            end
            reduction = rg.DimensionalityReduction.(level).(method).Reduction;
            metadata_object = rg.DimensionalityReduction.(level).ObjectGroup;
            [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, grouping_var);
            plot_cluster_outlines(rg,reduction, cluster_idx, nodeSz)
            if length(group_labels_comb) == 1
                group_labels_comb = strsplit(group_labels_comb, ' ');
            end
            legend(group_labels_comb,"Location","best","Box","off")
        end
        
        function plot_cluster_outlines(rg,reduction, cluster_idx, nodeSz, mapSz, sigma, cmap)
           arguments
               rg RecordingGroup
               reduction {isnumeric}
               cluster_idx (1,:) {isnumeric}
               nodeSz (1,1) {isnumeric} = 10
               mapSz {isnumeric} = 300
               sigma (1,1) {isnumeric} = mapSz/60
               cmap (:,3) {isnumeric} = othercolor('Spectral9',max(cluster_idx));
           end
           
           buffer_ratio = 0.5; 
           reduction = reduction - min(min(reduction)) * (1 + buffer_ratio);
           
           xrng = ...                      % range of x coordinates (with buffer)
               [min(reduction(:,1)) - buffer_ratio*range(reduction(:,1)),max(reduction(:,1)) + buffer_ratio*range(reduction(:,1))];
           yrng = ...                      % range of y coordinates (with buffer)
               [min(reduction(:,2)) - buffer_ratio*range(reduction(:,2)),max(reduction(:,2)) + buffer_ratio*range(reduction(:,2))];
           
           grid = zeros(mapSz,mapSz,max(cluster_idx));
           for i = 1:max(cluster_idx)                  % loop over all categories
               gridtemp = zeros(mapSz);        % temporary grid
               idx = cluster_idx == i;                 % get current category
               coori = reduction(idx,:);            % coordinates
               x = interp1(xrng,[1,mapSz],coori(:,1)); % map x coordinates to cells
               y = interp1(yrng,[1,mapSz],coori(:,2)); % same for y
               rx = round(x);                  % round to whole #
               ry = round(y);                  % ditto
               pts = (rx - 1)*mapSz + ry;
               gridtemp(pts) = 1;
               grid(:,:,i) = gridtemp;
           end
           reduction(:,1) = interp1(xrng,[1,mapSz],reduction(:,1));
           reduction(:,2) = interp1(yrng,[1,mapSz],reduction(:,2));
           
           sz = ceil(sqrt(mapSz));
           x = linspace(-sz,sz);
           y = linspace(-sz,sz);
           [mx,my] = meshgrid(x,y);
           kernel = exp(-(mx.^2 + my.^2)./(2*sigma.^2));
           kernel = kernel/sum(kernel(:));
           
           tot = zeros(mapSz,mapSz,3);
           colall = zeros(mapSz,mapSz,3,max(cluster_idx));
           for i = 1:max(cluster_idx)
               gridsmooth = conv2(grid(:,:,i),kernel,'same');
               prob = gridsmooth/max(gridsmooth(:));
               tot(:,:,i) = prob;
               
               P = bsxfun(@minus,ones(mapSz,mapSz,3),prob);
               C = ones(mapSz,mapSz,3);
               for k = 1:3
                   C(:,:,k) = cmap(i,k)*prob;
               end
               col = P + C;
               colall(:,:,:,i) = col;
           end
           
           [~,idx] = max(tot,[],3);
           cc = ones(mapSz,mapSz,3);
           for i = 1:max(cluster_idx)
               mask = idx == i;
               [u,v] = find(mask);
               for j = 1:length(u)
                   cc(u(j),v(j),:) = colall(u(j),v(j),:,i);
               end
           end
           
           for i = 1:3
               cc(:,:,i) = conv2(cc(:,:,i),kernel,'same');
           end
%            prct = round(mapSz*0.1);
%            cc(1:prct,:,:) = 1;
%            cc(:,1:prct,:) = 1;
%            cc(:,(end - prct + 1):end,:) = 1;
%            cc((end - prct + 1):end,:,:) = 1;
           ax = axes;
           hold(ax,'on');
           imagesc(cc)
           for i = 1:max(cluster_idx)
               idx = cluster_idx == i;
               coori = reduction(idx,:);scatter(coori(:,1),coori(:,2),nodeSz,repmat(cmap(i,:),sum(idx),1),'o','filled', 'MarkerEdgeColor','k');
           end
           xlim([mapSz*0.1 mapSz*0.9])
           ylim([mapSz*0.1 mapSz*0.9])
           axis off
%            axis image square;
        end
        
        function plot_single_cluster(rg, reduction, cluster_idx)
            arguments
                rg RecordingGroup
                reduction {isnumeric}
                cluster_idx (1,:) {isnumeric} = ones(1,size(reduction,1))
            end
            
            if length(cluster_idx)>100
                sz = 5;
            else
                sz = 20;
            end
            
            N_clust = length(unique(cluster_idx));
            
            nexttile
            if size(reduction,2) == 2
                scatter(reduction(:,1),reduction(:,2),sz,cluster_idx,'filled')
            else
                scatter3(reduction(:,1),reduction(:,2),reduction(:,3),sz,cluster_idx,'filled')
            end
            
            c_map = othercolor('Set19',N_clust);
            colormap(c_map)
            
        end
        
        function plot_dimensionality_reduction(rg,reduction,cluster_idx,group_labels_comb)
            arguments
                rg RecordingGroup
                reduction {isnumeric}
                cluster_idx (1,:) {isnumeric} = ones(1,size(reduction,1))
                group_labels_comb string = "Cluster " + unique(cluster_idx)
            end
            
            if length(cluster_idx)>500
                sz = 5;
            else
                sz = 20;
            end
            
            N_clust = length(unique(cluster_idx));
            figure('Color','w');
            if size(reduction,2) == 2
                nexttile
                scatter(reduction(:,1),reduction(:,2),sz,cluster_idx,'filled','MarkerEdgeColor','k')
                c_map = othercolor('Spectral9',N_clust);
                colormap(c_map)
                cb = colorbar; cb.Title.String = "Cluster Index";
                cb.Ticks = linspace(cb.Limits(1) + (N_clust-1)/(2*N_clust), cb.Limits(2) - (N_clust-1)/(2*N_clust), N_clust); cb.TickLabels = 1:N_clust;
                if N_clust > 1
                    xl = xlim; yl = ylim;
                    for i = 1:N_clust
                        nexttile
                        scatter(reduction(cluster_idx==i,1),reduction(cluster_idx==i,2),sz,c_map(i,:),'filled','MarkerEdgeColor','k')
                        title(group_labels_comb(i))
                        xlim(xl); ylim(yl)
                    end
                end
            else
                nexttile
                scatter3(reduction(:,1),reduction(:,2),reduction(:,3),sz,'filled','MarkerEdgeColor','k')
                c_map = othercolor('Spectral9',N_clust);
                colormap(c_map)
                cb = colorbar; cb.Title.String = "Cluster Index";
                cb.Ticks = linspace(cb.Limits(1) + (N_clust-1)/(2*N_clust), cb.Limits(2) - (N_clust-1)/(2*N_clust), N_clust); cb.TickLabels = 1:N_clust;
                if N_clust > 1
                    for i = 1:N_clust
                        nexttile
                        scatter3(reduction(cluster_idx==i,1),reduction(cluster_idx==i,2),reduction(cluster_idx==i,3),sz,c_map(i,:),'filled','MarkerEdgeColor','k')
                        title(group_labels_comb(i))
                    end
                end
            end
        end
        
        function plot_cluster_waveforms(rg,method)
            arguments
                rg RecordingGroup
                method string = "kmeans"
            end
            
            cluster_idx = [rg.Clustering.Unit.(method).Index];
            aligned_wf = rg.alignWaveforms([rg.Units]);
            N_clust = length(unique(cluster_idx));
            avg_wf = cell(1,N_clust);
            cluster_wfs = cell(1,N_clust);
            cluster_size = groupcounts(cluster_idx);
            time_vector = (1:size(aligned_wf,1))/20;
            
            for c = 1:N_clust
                cluster_wfs{c} = aligned_wf(:,cluster_idx == c);
                avg_wf{c} = median(cluster_wfs{c},2);
            end
            
            figure("Color","w");
            ax = nexttile;
            p = plot(time_vector,horzcat(avg_wf{:}),"LineWidth",2);
            ax.ColorOrder = othercolor('Spectral9',max(cluster_idx));
            
            l = legend("N = " + string(cluster_size)); l.Box = "off";
            box off; xlabel("Time [ms]"); ylim([-1 1])
            for c = 1:N_clust
                nexttile
                plot(time_vector,cluster_wfs{c},"LineWidth",0.1,'Color',[p(c).Color,0.1]);
                hold on
                plot(time_vector,horzcat(avg_wf{c}),"LineWidth",2,'Color',[0 0 0 0.5])
                box off; axis tight; xlabel("Time [ms]"); ylim([-1 1])
            end
        end
        
        function plot_feature_trajectories(rg, metadata_field, timepoints, feature_group, normalization, feature_name, grouping_var, useClustered, tolerance)
            arguments
                rg RecordingGroup
                metadata_field string
                timepoints (1,:) {isnumeric}
                feature_group string
                normalization string = "scaled" %"baseline" or "scaled"
                feature_name string = []
                grouping_var = "Mutation"
                useClustered logical = false
                tolerance = 1
            end
            switch feature_group
                case "UnitFeatures"
                    unit_features = "all";
                    network_features = [];
                case "NetworkFeatures"
                    unit_features = [];
                    network_features = "all";
            end
            [feature_table, culture_array] = aggregateCultureFeatureTables(rg, level, metadata_field, timepoints, tolerance, network_features, unit_features, normalization, useClustered);
            [group_idx, group_labels_comb] = combineMetadataIndices(rg, culture_array, grouping_var);
            N_features = size(feature_table,2)/numel(timepoints);
            N_groups = length(group_labels_comb);
            jitter = linspace(-1,1,N_groups);
            c = othercolor('RdBu4',N_groups);
            fontsz = 8;
            
            figure('Color','w','Position',[100 100 1500 1000]);
            tiledlayout('flow','TileSpacing','tight')
            for f = 1:N_features
                nexttile
                feature_parts = strsplit(feature_table.Properties.VariableNames{f},"_");
                feature_name = [feature_parts{1:end-1}];
                feature_matrix = feature_table(:,f:N_features:end).Variables;
                for g = 1:N_groups
                    x = timepoints + jitter(g);
                    xx = linspace(min(x),max(x),100);
                    data_mat = feature_matrix(group_idx==g,:);
                    y = mean(data_mat,1,'omitnan');
                    yy = pchip(x,y,xx);
                    y_err = std(data_mat,[],1,'omitnan');
                    plot(xx,yy,'Color',c(g,:),'HandleVisibility','off')
                    hold on
                    errorbar(timepoints+jitter(g),y,y_err,...
                        'LineWidth',1,'Color',c(g,:),'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',2,...
                        'MarkerFaceColor',c(g,:));
                    set(gca,'FontSize',fontsz)
                    marg = get(gca,'ylabel');
                    set(marg,'Margin',3)
                end
                title(feature_name)
                xlabel(metadata_field)
                box off
            end
            l = legend(group_labels_comb,'Box','off');
%             l.Position = [0.02 0.8 0.1 0.1];
        end
        
        function [value_array, group_labels_comb] = plot_cluster_densities(rg, level, method, grouping_var, n_bins)
            arguments
                rg RecordingGroup
                level string = "Unit" %Unit, Recording, Culture
                method string = "UMAP" %dim reduction method: "UMAP","PCA"
                grouping_var string = "Mutation"
                n_bins (1,1) {isnumeric} = 100
            end
            reduction = rg.DimensionalityReduction.(level).(method).Reduction;
            metadata_object = rg.DimensionalityReduction.(level).ObjectGroup;
            [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, grouping_var);
            N_clust = length(unique(cluster_idx));
            edges = {linspace(min(reduction(:,1)), max(reduction(:,1)),n_bins) linspace(min(reduction(:,2)), max(reduction(:,2)),n_bins)};
            value_array = zeros(n_bins,n_bins,N_clust);
            
            figure('Color','w');
            tiledlayout('flow','TileSpacing','compact')
            for i = 1:N_clust
                
                plot_idx = cluster_idx==i;
                data = [reduction(plot_idx,1),reduction(plot_idx,2)];
                ax(i) = nexttile;
                values = hist3(data, edges);
                value_array(:,:,i) = values;
%                 values = values.'./length(cluster_idx);
                imagesc(values);
                
                title(group_labels_comb(i))
                xticks([])
                yticks([])
            end
%             arrayfun(@(x) set(x,'CLim',[0 max(value_array,[],'all')]),ax)
        end
        
        function plot_cluster_shifts(rg, value_array, group_labels_comb)
            arguments
                rg RecordingGroup
                value_array double
                group_labels_comb string
            end
            N_clust = size(value_array,3);
            clust_comp = combvec(1:N_clust,1:N_clust);
            clust_comp = clust_comp(:,clust_comp(1,:) > clust_comp(2,:));
            N_data = sum(sum(value_array(:,:,1)));
            comp_array = zeros(size(value_array));
            
            for i = 1:size(clust_comp,2)
                comp_data = value_array(:,:,clust_comp(2,i)) - value_array(:,:,clust_comp(1,i));
                comp_array(:,:,i) = comp_data./N_data;
                ax(i) = nexttile;
                imagesc(comp_array(:,:,i))
                title(group_labels_comb(clust_comp(2,i)) + " -> " + group_labels_comb(clust_comp(1,i)))
                xticks([])
                yticks([])
            end
            max_change = max(abs(comp_array),[],'all')/5;
            cmap = othercolor('RdBu9',100);
            colormap(cmap)
            arrayfun(@(x) set(x,'CLim',[-max_change max_change]),ax)
        end
    end
end