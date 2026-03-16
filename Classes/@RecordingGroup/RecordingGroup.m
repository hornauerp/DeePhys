classdef RecordingGroup < handle
    
    properties
        Recordings
        Units
        Cultures
        Parameters % Inclusion/exclusion criteria for recording selection
        DimensionalityReduction
        Clustering
        Classification
        AppliedClassification
        Regression
    end
    
    methods (Hidden)
        function rg = RecordingGroup(recording_array, parameters)
            arguments
                recording_array MEArecording
                parameters struct = struct();
            end
            
            rg.parseParameters(parameters);
            rg.harmonizeMetadata(recording_array);
            keep_recordings = rg.filterRecordingArray([recording_array.Metadata]);
            rg.Recordings = recording_array(keep_recordings);
            rg.Units = [rg.Recordings.Units];
            if isfield([rg.Recordings.Metadata],'ChipID')
                rg.Cultures = rg.groupCultures();
            end
            fprintf('RecordingGroup: %i recordings, %i units, %i cultures\n', ...
            length(rg.Recordings), length(rg.Units), length(rg.Cultures));
        end
        
        function parseParameters(rg, parameters)
            rg.Parameters = rg.returnDefaultParams();
            if isempty(fieldnames(parameters))
                warning("No parameters provided, using default ones")
            else
                rg.Parameters = parseStructParameters(rg.Parameters, parameters);
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
        
        function harmonizeMetadata(rg, recording_array)
            md_fields = arrayfun(@(x) fields(x.Metadata),recording_array,'un',0);
            unique_fields = unique(string(vertcat(md_fields{:})));
            for r = 1:length(recording_array)
                missing_field = unique_fields(~matches(unique_fields,md_fields{r}));
                for mf = 1:length(missing_field)
                    recording_array(r).Metadata.(missing_field(mf)) = [];
                end
            end
        end
        
        function cultures = groupCultures(rg, grouping_var)
            arguments
               rg RecordingGroup
               grouping_var string = "" % Metadata field that varies across recordings within one culture
                                        % (e.g. "Concentration" for dose-response, "RecordingDate" for timepoints)
            end
            cultures = Culture.empty;
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
                   if grouping_var ~= ""
                       group_metadata = [culture_array.Metadata];
                       group_values = unique([group_metadata.(grouping_var)]);
                       for gv = group_values
                           inclusion = {{grouping_var, gv}};
                           group_idx = rg.filterRecordingArray(group_metadata, inclusion);
                           group_array = culture_array(group_idx);
                           cultures = [cultures, Culture(group_array, grouping_var, gv)]; %#ok<AGROW>
                       end
                   else
                       cultures = [cultures, Culture(culture_array)]; %#ok<AGROW>
                   end
                end
            end
        end
    end
    
    methods (Static)
       function defaultParams = returnDefaultParams()
           %Input to the filterRecordingArray function
           defaultParams.Selection.Inclusion = {}; % Cell array of cell arrays with {fieldname, value}; empty = include all recordings
           defaultParams.Selection.Exclusion = {}; % Cell array of cell arrays with {fieldname, value} to exclude
       end
       
       function aligned_wf = alignWaveforms(unit_array, target_sr)
       % ALIGNWAVEFORMS  Interpolate and trough-align reference waveforms.
       %   Delegates to FeatureAssembly.alignWaveforms for the actual logic.
           arguments
               unit_array Unit
               target_sr  (1,1) double = 0
           end
           aligned_wf = FeatureAssembly.alignWaveforms(unit_array, target_sr);
       end
       
       function [clf, train_acc] = create_classifier(X_train, Y_train, alg, N_hyper)
       % CREATE_CLASSIFIER  Train a classifier. Delegates to MLPipeline.createClassifier.
           [clf, train_acc] = MLPipeline.createClassifier(X_train, Y_train, alg, N_hyper);
       end
       
       function [Y_train, Y_test, train_idx, test_idx] = cv_split(Y, cv, k)
       % CV_SPLIT  Extract train/test split for fold k. Delegates to MLPipeline.cvSplit.
           [Y_train, Y_test, train_idx, test_idx] = MLPipeline.cvSplit(Y, cv, k);
       end

       function [new_group_idx, new_group_labels] = poolMetadataValues(group_idx, group_labels, classification_val)
       % POOLMETADATAVALUES  Pool metadata values. Delegates to MLPipeline.poolMetadataValues.
           [new_group_idx, new_group_labels] = MLPipeline.poolMetadataValues(group_idx, group_labels, classification_val);
       end
       
       function [reduction, cmap] = plot_cluster_outlines(reduction, cluster_idx, ax, plot_centroid, nodeSz, mapSz, sigma, cmap)
           arguments
               reduction {isnumeric}
               cluster_idx (1,:) {isnumeric}
               ax = axes;
               plot_centroid = false
               nodeSz (1,1) {isnumeric} = 10
               mapSz {isnumeric} = 500
               sigma (1,1) {isnumeric} = mapSz/6
               cmap (:,3) {isnumeric} = othercolor('Spectral6',max(cluster_idx)); %othercolor("Set29",4);%
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
           
           hold(ax,'on');
           imagesc(cc)
           for i = 1:max(cluster_idx)
               idx = cluster_idx == i;
               coori = reduction(idx,:);
               if plot_centroid
                   s = scatter(coori(:,1),coori(:,2),nodeSz,repmat(cmap(i,:),sum(idx),1),'o','filled', 'MarkerFaceAlpha',0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.25);
                   s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                   scatter(mean(coori(:,1)),mean(coori(:,2)),nodeSz*2,cmap(i,:),'o','filled', 'MarkerEdgeColor','k')
               else
                   scatter(coori(:,1),coori(:,2),nodeSz,repmat(cmap(i,:),sum(idx),1),'o','filled', 'MarkerEdgeColor','k');
               end
           end
           xlim([mapSz*0.1 mapSz*0.9])
           ylim([mapSz*0.1 mapSz*0.9])
           axis off
           %            axis image square;
       end
       
    end
    
    methods 
        
        
        function runANOVA(rg, level, metadata_field, values, tolerance, network_features, unit_features, normalization, useClustered)
           arguments
               rg RecordingGroup
               level string = "Recording"
               metadata_field string = "DIV"
               values {isnumeric} = [7,14,21,28] %refers to values of metadata_field; NaN to use the maximum number of available timepoints
               tolerance {isnumeric} = 1 %deviation from the actual age that will still be considered valid (e.g. age = 7 and tolerance = 1 considers DIVs 6-8)
               network_features string = "all"
               unit_features string = "all"
               normalization string = "baseline" % "baseline" (divided by first timepoint) or "scaled (scaled between 0 and 1)
               useClustered logical = false
           end
           
           [feature_table, culture_array, values] = aggregateCultureFeatureTables(rg, level, metadata_field, values, tolerance, network_features, unit_features, normalization, useClustered);
           %%% TO DO %%%
        end
        
        
        
        
        
        
        
                
        function [metrics_special, train_acc, avg_pred_imp, conf_mat] = assessClassifier(rg, clf_result)
            arguments
               rg RecordingGroup
               clf_result struct
            end
            try
                y = vertcat(clf_result.Y_test);
            catch
                y = horzcat(clf_result.Y_test)';
            end
            y_hat = vertcat(clf_result.Y_pred);
            
            conf_mat = confusionmat(y,y_hat);
            metrics_special = multiclass_metrics_special(conf_mat);

            pred_imp = vertcat(clf_result.predImp);
            
            
            train_acc = mean([clf_result.train_acc]);
            avg_pred_imp = mean(pred_imp);
            
            fprintf('Training accuracy: %.2f\n', train_acc)
            fprintf('F1-score: %.2f\n', metrics_special.F1_score)
        end
        
        
        function [train_acc, test_acc, avg_score] = classifyByFeatureGroupsAndGroupingVar(rg, level, alg, classification_var, pooling_vals, useClustered,...
                grouping_var, grouping_value_list, normalization, normalization_var, N_hyper, K_fold, tolerance)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                alg string = "rf" %rf, svm,
                classification_var string = "Mutation" %Metadata field that determines y_test/y_train
                pooling_vals cell = {} %Values corresponding to classification_var that is being classified for(e.g. "LNA")
                useClustered logical = false
                grouping_var string = "DIV" %Metadata field that groups recordings to cultures
                grouping_value_list = nan %Selected values corresponding to grouping_var
                normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
                normalization_var string = "PlatingDate" % normalization by each value of normalization_var
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
                tolerance double = 1
            end
            % Find all existing feature names
             all_unit_features = {"ActivityFeatures", "WaveformFeatures", ["RegularityFeatures","Catch22"], "all"};%[string(fieldnames([rg.Recordings.NetworkFeatures])); "all"];
             all_network_features = {["Regularity", "Catch22"], "Burst", "GraphFeatures", "all"};%[string(fieldnames([rg.Recordings.UnitFeatures])); "all"];
            
            train_acc = nan(length(all_network_features) + length(all_unit_features) + 1, length(grouping_value_list));
            test_acc = nan(length(all_network_features) + length(all_unit_features) + 1, length(grouping_value_list));
            avg_score = nan(length(all_network_features) + length(all_unit_features) + 1, length(grouping_value_list));
            
            % First iterate over network features
            unit_features = [];
            for iNW = 1:length(all_network_features)
                for iGV = 1:length(grouping_value_list)
                    result = classifyByFeatureGroups(rg, level, alg, classification_var, pooling_vals, all_network_features{iNW}, unit_features, [], useClustered,...
                        grouping_var, grouping_value_list(iGV), normalization, normalization_var, N_hyper, K_fold, tolerance);
                    [metrics, train_acc(iNW, iGV)] = assessClassifier(rg, result);
                    test_acc(iNW, iGV) = metrics.F1_score;
                    avg_score(iNW, iGV) = metrics.Accuracy;
                end
            end
            
            % Then over all unit features
            network_features = [];
            for iU = 1:length(all_unit_features)
                for iGV = 1:length(grouping_value_list)
                    result = classifyByFeatureGroups(rg, level, alg, classification_var, pooling_vals, network_features, all_unit_features{iU}, [], useClustered,...
                        grouping_var, grouping_value_list(iGV), normalization, normalization_var, N_hyper, K_fold, tolerance);
                    [metrics, train_acc(iU + iNW, iGV)] = assessClassifier(rg, result);
                    test_acc(iU + iNW, iGV) = metrics.F1_score;
                    avg_score(iU + iNW, iGV) = metrics.Accuracy;
                end
            end
            
            % Finally over all feature groups combined
            network_features = "all";
            unit_features = "all";
            for iGV = 1:length(grouping_value_list)
                result = classifyByFeatureGroups(rg, level, alg, classification_var, pooling_vals, network_features, unit_features, [], useClustered,...
                    grouping_var, grouping_value_list(iGV), normalization, normalization_var, N_hyper, K_fold, tolerance);
                [metrics, train_acc(end, iGV)] = assessClassifier(rg, result);
                test_acc(end, iGV) = metrics.F1_score;
                avg_score(end, iGV) = metrics.Accuracy;
            end
            rg.Classification.FeatureGroupsByGroupingVar.TrainAcc = train_acc;
            rg.Classification.FeatureGroupsByGroupingVar.TestAcc = test_acc;
            rg.Classification.FeatureGroupsByGroupingVar.AvgScore = avg_score;
            rg.Classification.FeatureGroupsByGroupingVar.FeatureGroups = [all_network_features all_unit_features {"combined"}];
            rg.Classification.FeatureGroupsByGroupingVar.GroupingVar = grouping_var;
            rg.Classification.FeatureGroupsByGroupingVar.GroupingValueList = grouping_value_list;
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
        
        function [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, metadata_names, pooling_vals)
            arguments
                rg RecordingGroup
                metadata_object %Can be a string referring to the respective group (Unit,Recording,Culture) within the RecordingGroup, or a Unit,Recording or Culture array
                metadata_names string
                pooling_vals cell = {}
            end
            
            idx_array = cell(1,length(metadata_names));
            group_array = cell(1,length(metadata_names));
            for m = 1:length(metadata_names)
                [idx, group] = returnMetadataArray(rg, metadata_object, metadata_names(m)); %{end-m+1}
                if ~isempty(pooling_vals)
                    new_idx_array = [];
                    new_groups = string();
                    for p = 1:length(pooling_vals{m})
                        pool_idx = find(matches(string(group), string(pooling_vals{m}{p})));
                        new_idx_array(p,:) = matches(string(idx), string(pool_idx)) * p;
                        new_groups(p) = join(string(group(pool_idx)),' ');
                    end
                    idx_array{end-m+1} = sum(new_idx_array);
                    group_array{end-m+1} = new_groups;
                else
                    idx_array{end-m+1} = idx;
                    group_array{end-m+1} = group;
                end
            end
            group_vals = cellfun(@unique, idx_array, "un",0);
            group_combs = combvec(group_vals{:});
            group_idx = vertcat(idx_array{:});
            [~,cluster_idx] = ismember(group_idx', group_combs','rows');
            group_labels = arrayfun(@(x) group_array{x}(group_combs(x,:))',1:length(group_array),"un",0);
            if length(metadata_names) == 1
                group_labels_comb = [group_labels{:}];
            else
                group_labels_comb = join(string([group_labels{end:-1:1}]));
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
        
        function culture = findCulture(rg, rec)
            arguments
                rg RecordingGroup
                rec MEArecording
            end
            culture_idx = cellfun(@(x) any(ismember(x, rec)),rg.Cultures);
            culture = rg.Cultures{culture_idx};
        end

        function [sel_recs, group_labels_comb] = findBestRecordings(rg, metadata_names, pooling_vals)
            arguments
                rg RecordingGroup
                metadata_names string
                pooling_vals cell = {}
            end
            if ~isempty(metadata_names)
                [cluster_idx, group_labels_comb] = rg.combineMetadataIndices(rg.Recordings, metadata_names, pooling_vals);
            else
                cluster_idx = ones(1,length(rg.Recordings));
                group_labels_comb = [];
            end
            N_units = cellfun(@length, {rg.Recordings.Units});

            sel_recs = MEArecording();
            for c = 1:length(unique(cluster_idx))
                sel_idx = find(cluster_idx == c);
                max_idx = sel_idx(N_units(sel_idx) == max(N_units(sel_idx)));
                sel_recs(c) = rg.Recordings(max_idx(1));
            end
        end

        function feature_names = returnFeatureNames(rg, feature_groups, feature_subgroups)
            arguments
               rg RecordingGroup
               feature_groups string = ["UnitFeatures","NetworkFeatures"]
               feature_subgroups cell  = {}
            end
            
            
            for fg = 1:length(feature_groups)
                if isempty(feature_subgroups)
                    feature_sg = string(fieldnames([rg.Recordings.(feature_groups(fg))]));
                else
                    feature_sg = feature_subgroups{fg};
                end
                features = [];
                for fs = 1:length(feature_sg)
                    features = [features string(rg.Recordings(1).(feature_groups(fg)).(feature_sg(fs)).Properties.VariableNames)];
                end
                feature_names{fg} = features;
            end
        end
    
        function [cluster_idx, cmap] = clusterByFeatures(rg, input_type, level, method, N_clust, grouping_var)
           arguments
                rg RecordingGroup
                input_type string = "UMAP" %"RAW","UMAP","PCA"
                level string = "Unit" %Unit, Recording, Culture
                method string = "kmeans" %kmeans, hierarchical, spectral, gmm, louvain
                N_clust {isnumeric} = nan %Number of imposed clusters; nan selects N_clust by finding the optimal silhouette score
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
                   
               case "louvain"
                   assert(input_type == "UMAP","UMAP graph needed for louvain clustering")
                   [M,Q] = community_louvain(rg.DimensionalityReduction.(level).UMAP.Graph, N_clust); %N_Clust refers here to the resolution parameter
                   
           end
           if isnan(N_clust)
               evaluation = evalclusters(feat_mat,cluster_fun,'silhouette',"KList",1:10);
               N_clust = evaluation.OptimalK;
               cluster_idx = cluster_fun(feat_mat,N_clust);
           elseif method == "louvain"
               cluster_idx = M;
%                fprintf('Modularity: %.2f\n',Q)
           else
               cluster_idx = cluster_fun(feat_mat,N_clust);
           end 
           rg.Clustering.(level).(method).Index = cluster_idx;
           rg.Clustering.(level).(method).k = N_clust;
           rg.Clustering.(level).(method).Input = input_type;
           
%            figure('Color','w');
           [reduction, cmap] = RecordingGroup.plot_cluster_outlines(feat_mat, cluster_idx, gca);
           title(sprintf('%s clustering on %ss',method, level))
        end
        
        function cluster_purity = calculateClusterPurity(rg, input_type, method, cluster_vars, pooling_vals)
           arguments
               rg RecordingGroup
               input_type string = "UMAP"
               method string = "kmeans" %Clustering method
               cluster_vars string = "Mutation" %Metadata variable for the ground truth
               pooling_vals cell = {}
           end
           
           [ground_truth_labels, ~] = combineMetadataIndices(rg, "Culture", cluster_vars, pooling_vals);
           N_clust = max(ground_truth_labels);
           if method == "louvain"
               predicted_labels = 0;
               res_param = 0.9;
               while max(predicted_labels) < N_clust
                   predicted_labels = rg.clusterByFeatures(input_type, "Recording", method, res_param);
                   res_param = res_param + 0.01;
                   if max(predicted_labels) > N_clust
                       fprintf('Resolution parameter: %.2f ',res_param)
                       error('Could not find suitable resolution parameter')
                   end
               end
           else
               predicted_labels = rg.clusterByFeatures(input_type, "Recording", method, N_clust);
           end
           cm = confusionmat(ground_truth_labels,predicted_labels);

           cluster_purity = sum(max(cm)) / length(predicted_labels);
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
            rg.Units = [rg.Recordings.Units];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [cluster_idx, group_labels_comb] = plot_true_clusters(rg, level, method, grouping_var, options)
            arguments
                rg RecordingGroup
                level string = "Unit" %Unit, Recording, Culture
                method string = "UMAP" %dim reduction method: "UMAP","PCA"
                grouping_var string = "Mutation"
                options.pooling_vals cell = {}
                options.plot_centroid = false
                options.nodeSz (1,1) {isnumeric} = 10
                options.mapSz {isnumeric} = 300
                options.sigma (1,1) {isnumeric} = 5
                options.cmap (:,3) double = []
            end
            reduction = rg.DimensionalityReduction.(level).(method).Reduction;
            metadata_object = rg.DimensionalityReduction.(level).ObjectGroup;
            [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, grouping_var, options.pooling_vals);
            if isempty(options.cmap)
                options.cmap = othercolor('Mrainbow',max(cluster_idx));%Mdarkrainbow
            end
            ax = gca;
            RecordingGroup.plot_cluster_outlines(reduction, cluster_idx, ax, options.plot_centroid, options.nodeSz, options.mapSz, options.sigma, options.cmap);
            if length(group_labels_comb) == 1
                group_labels_comb = strsplit(group_labels_comb, ' ');
            end
            legend(group_labels_comb,"Location","best","Box","off")
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
            c_map = othercolor('Spectral9',N_clust);
            figure('Color','w');
            if size(reduction,2) == 2
                nexttile
                scatter(reduction(:,1),reduction(:,2),sz,cluster_idx,'filled','MarkerEdgeColor','k')
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
                scatter3(reduction(:,1),reduction(:,2),reduction(:,3),sz,cluster_idx,'filled','MarkerEdgeColor','k')
                colormap(c_map)
                cb = colorbar; cb.Title.String = "Cluster Index";
                cb.Ticks = linspace(cb.Limits(1) + (N_clust-1)/(2*N_clust), cb.Limits(2) - (N_clust-1)/(2*N_clust), N_clust); cb.TickLabels = 1:N_clust;
                if N_clust > 1
                    xl = xlim; yl = ylim; zl = zlim;
                    for i = 1:N_clust
                        nexttile
                        scatter3(reduction(cluster_idx==i,1),reduction(cluster_idx==i,2),reduction(cluster_idx==i,3),sz,c_map(i,:),'filled','MarkerEdgeColor','k')
                        title(group_labels_comb(i))
                        xlim(xl); ylim(yl);zlim(zl)
                    end
                end
            end
        end
        
        function [colors, avg_wf] = plot_cluster_waveforms(rg,method,cluster_idx, object_group)
            arguments
                rg RecordingGroup
                method string = "kmeans"
                cluster_idx {isnumeric} = []
                object_group = []
            end
            if isempty(cluster_idx)
                cluster_idx = [rg.Clustering.Unit.(method).Index];
            end

            if isempty(object_group)
                aligned_wf = rg.alignWaveforms([rg.Units]);
            else
                aligned_wf = rg.alignWaveforms(object_group);
            end
            N_clust = length(unique(cluster_idx));
            avg_wf = cell(1,N_clust);
            cluster_wfs = cell(1,N_clust);
            cluster_size = histcounts(cluster_idx);
            time_vector = (1:size(aligned_wf,1))/20;
            
            for c = 1:N_clust
                cluster_wfs{c} = aligned_wf(:,cluster_idx == c);
                avg_wf{c} = median(cluster_wfs{c},2);
            end

            avg_time_vector = linspace(min(time_vector),max(time_vector),max([200,length(time_vector)]));
            smoothed_avg = pchip(time_vector,horzcat(avg_wf{:})',avg_time_vector);

            figure("Color","w");
            ax = nexttile;

            p = plot(avg_time_vector,smoothed_avg,"LineWidth",2);
            ax.ColorOrder = othercolor('Spectral9',max(cluster_idx));
            
            l = legend("N = " + string(cluster_size)); l.Box = "off";
            box off; xlabel("Time [ms]"); ylim([-1 1])
            for c = 1:N_clust
                nexttile
                plot(time_vector,cluster_wfs{c},"LineWidth",0.1,'Color',[p(c).Color,0.1]);
                hold on
                plot(avg_time_vector,smoothed_avg(c,:),"LineWidth",2,'Color',[0 0 0 0.5])
                box off; axis tight; xlabel("Time [ms]"); ylim([-1 1])
            end
            colors = vertcat(p.Color);
        end
        
        function [feature_table, mean_mat, sd_mat, group_labels_comb] = plot_feature_trajectories(rg, level, grouping_var, grouping_values, network_features, unit_features, normalization,...
                feature_names, comp_var, pooling_vals, useClustered, tolerance, colors)
            arguments
                rg RecordingGroup
                level string
                grouping_var string
                grouping_values (1,:) {isnumeric}
                network_features string
                unit_features string
                normalization string = [] %"baseline" or "scaled"
                feature_names string = []
                comp_var string = "Mutation"
                pooling_vals cell = {}
                useClustered logical = false
                tolerance = 0
                colors double =[]
            end
            [feature_table, culture_array, grouping_values] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            
            if ~isempty(feature_names)
                sel_idx = startsWith(feature_table.Properties.VariableNames, feature_names);
                if isempty(sel_idx)
                    error("Could not find " + feature_names)
                end
                feature_table = feature_table(:,sel_idx);
            end
            
            grouping_values = 1:length(grouping_values);
            if normalization == "baseline"
               % grouping_values =  grouping_values(2:end);
               grouping_values = 1:(length(grouping_values)-1);
            end
            
            N_features = size(feature_table,2)/numel(grouping_values);
            
            if level == "Recording"
                [group_idx, group_labels_comb] = combineMetadataIndices(rg, culture_array, comp_var, pooling_vals);
                N_groups = length(group_labels_comb);
            elseif level == "Unit"
                %                 group_idx = [rg.Units.ClusterID];
                group_idx = [rg.DimensionalityReduction.Unit.ObjectGroup.ClusterID];
                N_groups = length(unique(group_idx));
                group_labels_comb = "Cluster " + [1:N_groups];
            else
                error("Unknown level")
            end
            
            fontsz = 7;
            
            if isempty(colors)
               colors = othercolor('RdBu4',N_groups);
            end
            
            min_diff = min([1, diff(grouping_values)]); % In case we only plot one value
            jitter = linspace(-min_diff/10,min_diff/10,N_groups);%jitter = zeros(1,N_groups); %Change to have jitter
            
            mean_mat = nan(N_features, length(grouping_values), N_groups);
            sd_mat = nan(N_features, length(grouping_values), N_groups);
            features = string();
            
            figure('Color','w','Position',[100 100 1500 1000]);
            tiledlayout('flow','TileSpacing','tight')
            for f = 1:N_features
                if ~isempty(feature_names)
                    sel_idx = find(startsWith(feature_table.Properties.VariableNames, feature_names(f)));
                else
                    sel_idx = f:N_features:size(feature_table,2);
                end
                feature_parts = strsplit(feature_table.Properties.VariableNames{sel_idx(1)},"_");
                features(f) = [feature_parts{1:end-1}];
                feature_matrix = feature_table(:,sel_idx).Variables;
                
                nexttile
                
                for g = 1:N_groups
                    x = grouping_values + jitter(g);
                    
                    data_mat = feature_matrix(group_idx==g,:);
                    data_mat(isoutlier(data_mat,'ThresholdFactor',10)) = nan;

                    
                    if normalization == "baseline"
                        x = (1:(length(grouping_values)+1)) + jitter(g);
                        plot_x = 1:length(x);
                        log_data = log(data_mat);
                        log_data(isinf(log_data)) = nan;
                        y = [0 mean(log_data,1,'omitnan')];
                        y_err = [0 std(log_data,[],1,'omitnan')];
                    else
                        y = mean(data_mat,1,'omitnan');
                        y_err = std(data_mat,[],1,'omitnan');
                        plot_x = grouping_values;
                    end
                    % if length(x) > 1
                    %     xx = linspace(min(x),max(x),((max(x) - min(x))*5)+1);
                    %     yy = makima(x,y,xx);
                    % else
                    % 
                    % end
                    
                    plot(x,y,'Color',colors(g,:),'HandleVisibility','off')
                    hold on
                    errorbar(plot_x+jitter(g),y,y_err,...
                        'LineWidth',1,'Color',colors(g,:),'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',5,...
                        'MarkerFaceColor',colors(g,:));
                    set(gca,'FontSize',fontsz)
                    marg = get(gca,'ylabel');
                    set(marg,'Margin',3)
                    if normalization == "baseline"
                        mean_mat(f,:,g) = y(2:end);
                        sd_mat(f,:,g) = y_err(2:end);
                    else
                        mean_mat(f,:,g) = y;
                        sd_mat(f,:,g) = y_err;
                    end
                end
                l = legend(group_labels_comb,'Box','off');
                
                title(features(f))
                xlabel(grouping_var)
                box off
            end
            
            %             l.Position = [0.02 0.8 0.1 0.1];
        end
        
        function [color_mat, features] = plot_feature_heatmap(rg, level, grouping_var, grouping_values, network_features, unit_features, feature_names, normalization, comp_var, pooling_vals, useClustered, tolerance, color_lim)
           arguments
                rg RecordingGroup
                level string
                grouping_var string
                grouping_values (1,:) {isnumeric}
                network_features string
                unit_features string
                feature_names string
                normalization string = [] %[], "baseline" or "scaled"
                comp_var = "Mutation"
                pooling_vals cell = {}
                useClustered logical = false
                tolerance = 1
                color_lim double = 3 %Maximum value to cap colormap
           end
           
           [feature_table, culture_array, grouping_values] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
           
           if ~isempty(feature_names)
                vars = string(feature_table.Properties.VariableNames);
                sel_vars = startsWith(vars, feature_names);
                feature_table = feature_table(:,sel_vars);
           end
            
           %We return grouping_values in case we extracted the values by
           %setting it to nan 
           [group_idx, group_labels_comb] = combineMetadataIndices(rg, culture_array, comp_var, pooling_vals);
           features = [];
           feat_mat = feature_table.Variables;
           feat_mat(isinf(feat_mat)) = nan;
           norm_mat = normalize(feat_mat,'range',[0 1]);
           means = arrayfun(@(x) mean(norm_mat(group_idx == x,:),'omitnan'),1:length(unique(group_idx)),'un',0);
           color_vector = log(means{1}./means{2});
           color_mat = reshape(color_vector,[],length(grouping_values));
           for f = 1:length(feature_table.Properties.VariableNames)
               feature_parts = strsplit(feature_table.Properties.VariableNames{f},"_");
               feat = string(strjoin(feature_parts(1:end-1),'_'));
               if isempty(features) || ~matches(feat, features)
                   features = [features feat];
               end
           end
           
           if ~isempty(feature_names) %Sort features to match input
               [~,idx] = ismember(feature_names, features);
               color_mat = color_mat(idx,:);
               features = feature_names;
           end
           
           if normalization == "baseline"
              color_mat(:,1) = [];
              grouping_values(1) = [];
           end
           
           figure('Color','w');
           imagesc(color_mat)
           cm = othercolor('RdBu9',100);
           colormap(cm([1:45,55:end],:)) %Make smaller differences more visible
           clim = [log(1/color_lim) log(color_lim)];
           set(gca,'CLim',clim) %Set color limits to 1/3 and 3 
           xticks(1:size(color_mat,2))
           xticklabels(grouping_values)
           xlabel(grouping_var)
           yticks(1:size(color_mat,1))
           yticklabels(features)
           cb = colorbar;
           cb.Location = 'northoutside';
           cb.Title.String = group_labels_comb(1) + "/" + group_labels_comb(2);
           cb.Title.FontWeight = 'bold';
           cb.Ticks = [clim(1) 0 clim(end)];
           cb.TickLabels = [sprintf("<%.1f",1/color_lim) + color_lim, "1", ">" + color_lim];
           set(gca,'FontSize',7)
           set(gca,'TickLength',[0 0])
        end
        
        function [value_array, group_labels_comb] = plot_cluster_densities(rg, level, method, grouping_var, pooling_vals, n_bins, smoothing_factor)
            arguments
                rg RecordingGroup
                level string = "Unit" %Unit, Recording, Culture
                method string = "UMAP" %dim reduction method: "UMAP","PCA"
                grouping_var = "Mutation" %Metadata name or cluster index
                pooling_vals cell = {}
                n_bins (1,1) {isnumeric} = 50
                smoothing_factor double = 0.5
            end
            reduction = rg.DimensionalityReduction.(level).(method).Reduction;
            if isstring(grouping_var)
                metadata_object = rg.DimensionalityReduction.(level).ObjectGroup;
                [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, grouping_var, pooling_vals);
            else
                cluster_idx = grouping_var;
                group_labels_comb = string(unique(cluster_idx));
            end
            N_clust = length(unique(cluster_idx));
            edges = {linspace(min(reduction(:,1)), max(reduction(:,1)),n_bins) linspace(min(reduction(:,2)), max(reduction(:,2)),n_bins)};
            value_array = zeros(n_bins,n_bins,N_clust);
            [X,Y] = meshgrid(edges{1},edges{2});
            figure('Color','w');
            tiledlayout('flow','TileSpacing','compact')
            for i = 1:N_clust
                
                plot_idx = cluster_idx==i;
                data = [reduction(plot_idx,1),reduction(plot_idx,2)];
                ax(i) = nexttile;
                values = hist3(data, edges);
                embedded_values = zeros(size(values)*1.5);
                embed_coor = round(n_bins*0.25)+1:round(n_bins*1.25);
                embedded_values(embed_coor, embed_coor) = values;
                value_array(:,:,i) = values;
                % values = values.'./length(cluster_idx);
% surf(X,Y,values)
                smoothed_data = smoothdata2(embedded_values,"gaussian",SmoothingFactor=smoothing_factor);
                contourf(smoothed_data,'LineWidth',0.25,'FaceAlpha',0.5);
                axis square
                title(group_labels_comb(i))
                xticks([])
                yticks([])
                colormap(flipud(hot))
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
                comp_data = normalize(value_array(:,:,clust_comp(2,i)),'range',[0 1]) - normalize(value_array(:,:,clust_comp(1,i)),'range',[0 1]);
                comp_array(:,:,i) = comp_data./N_data;
                ax(i) = nexttile;
                imagesc(comp_array(:,:,i))
                title(group_labels_comb(clust_comp(2,i)) + " -> " + group_labels_comb(clust_comp(1,i)))
                xticks([])
                yticks([])
            end
            max_change = max(abs(comp_array),[],'all');
            cmap = othercolor('RdBu9',100);
            colormap(cmap)
            arrayfun(@(x) set(x,'CLim',[-max_change max_change]),ax)
        end
        
        function plot_regression_results(rg, regression_var, color_var, color_order)
            arguments
                rg RecordingGroup
                regression_var string
                color_var string = "Mutation"
                color_order double = []
            end
            y_pred = vertcat(rg.Regression.(regression_var).Y_pred);
            y_test = vertcat(rg.Regression.(regression_var).Y_test);
            reg_objects = [rg.Regression.(regression_var).objects];
            [color_idx, color_labels] = rg.returnMetadataArray(reg_objects, color_var);
            
            min_diff = min(diff(unique(y_test)));
            x_vals = y_test/min_diff;
            
            figure('Color','w');
            b = boxchart(x_vals, y_pred,'GroupByColor',color_idx,'MarkerSize',1);
            hold on
            p = plot([0 max(x_vals)*1.2],[0 max(y_test)*1.2],'k--'); p.Color(4) = 0.3;
            xticks(unique(x_vals))
            xticklabels(unique(y_test))
            xlabel(regression_var)
            ylim([0 max(y_test) + 0.5*min(y_test)])
            yticks(unique(y_test))
            ylabel(" Predicted " + regression_var)
            set(gca,'FontSize',8)
            leg = legend(color_labels);
            leg.Box = 'off';
            leg.Location = 'best';
            if isempty(color_order)
               color_order = othercolor('RdBu9',length(b)); 
            end
            arrayfun(@(x) set(b(x),'BoxFaceColor', color_order(x,:)),1:length(b))
        end
        
        function rel_counts = plot_cluster_proportions(rg, separation_var, pooling_vals, clust_method, colors)
            arguments
                rg RecordingGroup
                separation_var string
                pooling_vals cell = {}
                clust_method string = "spectral"
                colors double = []
            end
            [idx, labels] = rg.combineMetadataIndices(rg.DimensionalityReduction.Unit.ObjectGroup,separation_var, pooling_vals);
            cluster_idx = [rg.Clustering.Unit.(clust_method).Index];
            
            cluster_counts_mat = nan(length(labels),length(unique(cluster_idx)));
            for i = 1:length(labels)
                cluster_counts_mat(i,:) = histcounts(cluster_idx(idx == i),'BinLimits',[1 max(cluster_idx)]);
            end
            
            rel_counts = cluster_counts_mat./sum(cluster_counts_mat,2);

            b = bar(rel_counts,'stacked','FaceColor','flat');
            if ~isempty(colors)
                arrayfun(@(x) set(b(x),'CData',repmat(colors(x,:),[length(labels) 1])),1:length(b))
            end
            xlabel(separation_var)
            xticklabels(labels)
            ylabel("Cluster percentage")
            box off
            axis tight
        end
        
        % plot_group_feature_heatmap — placeholder removed, implement when needed

        
        function plot_unit_cluster_features(rg, cluster_idx, feature_names, colors)
            % Plots box plots per cluster for the selected features
            arguments
                rg RecordingGroup
                cluster_idx double
                feature_names string
                colors double = []
            end
            input_table = rg.Recordings.getUnitFeatures(rg.Units,"all");
            figure('Color','w');
            tiledlayout('flow')
            for f = 1:length(feature_names)
                x = input_table(:,feature_names(f)).Variables;
                nexttile
                b = boxchart(x,'GroupByColor',cluster_idx,'MarkerSize',1);
                title(feature_names(f))
                if ~isempty(colors)
                    arrayfun(@(x) set(b(x), 'BoxFaceColor',colors(x,:),'MarkerColor',colors(x,:)),1:length(b))
                end
                ylim([quantile(x,0.05),quantile(x,0.95)])
            end
            xticklabels([])
        end
        
        function color_mat = plot_unit_cluster_heatmap(rg, cluster_idx, feature_names, cmap)
            % Plots box plots per cluster for the selected features
            arguments
                rg RecordingGroup
                cluster_idx double
                feature_names string = []
                cmap double = othercolor('RdBu9',100);
            end
            N_Clust = length(unique(cluster_idx));
            
            input_table = rg.Recordings.getUnitFeatures(rg.Units,"all");
            if isempty(feature_names)
                feature_names = string([input_table.Properties.VariableNames]);
            end
            color_mat = nan(length(feature_names), N_Clust);
            feat_mat = normalize(input_table(:,feature_names).Variables);
            for c = 1:N_Clust
                color_mat(:,c) = mean(feat_mat(cluster_idx == c,:),'omitnan');
            end
            figure('Color','w');
            imagesc(color_mat)
            colormap(cmap)
            xticks(1:N_Clust)
            xticklabels(1:N_Clust)
            xlabel('Cluster ID')
            yticks(1:length(feature_names))
            yticklabels(feature_names)
            colorbar
            set(gca,'CLim',[-max(max(abs(color_mat))) max(max(abs(color_mat)))])
            set(gca,"TickLength",[0,0])
            setAllFontSizes(gcf,7)
        end
    end
end