classdef RecordingGroup < handle
    
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
                for pf = parameter_fields
                    if ismember(pf,fieldnames(rg.Parameters))
                        subfields = fieldnames(parameters.(pf{:}));
                        for sf = subfields
                            if ismember(sf,fieldnames(rg.Parameters.(pf{:})))
                                rg.Parameters.(pf{:}).(sf{:}) = parameters.(pf{:}).(sf{:});
                            else
                                warning("Unrecognized parameter field: %s.%s",pf{:},sf{:})
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
                    value_vector = {metadata_array.(inclusion{i}{1})};
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
    end
    
    methods 
        function feature_table = aggregateRecordingFeatureTables(rg, unit_features, network_features)
            arguments
                rg RecordingGroup
                unit_features string = "all"
                network_features string = "all"
            end
            
            feature_array = cell(1,length(rg.Recordings));
            for r = 1:length(rg.Recordings)
                feature_array{r} = getFeatureTable(rg.Recordings(r), unit_features, network_features);
            end
            feature_table = vertcat(feature_array{:});
        end
        
        function [feature_table, culture_array] = aggregateCultureFeatureTables(rg, age, tolerance, unit_features, network_features)
            arguments
                rg RecordingGroup
                age {isnumeric} = [7,14,21,28] %refers to DIV // 0 to use the maximum number of available timepoints
                tolerance {isnumeric} = 1 %deviation from the actual age that will still be considered valid (e.g. age = 7 and tolerance = 1 considers DIVs 6-8)
                unit_features string = "all"
                network_features string = "all"
            end
            
            rec_metadata = [rg.Recordings.Metadata];
            ages = unique([rec_metadata.DIV]);
            age_dev = abs(age - ages');
            age_count = sum(age_dev<=tolerance);
            if any(age_count == 0)
               warning("DIV " + num2str(age(age_count == 0)) + " without corresponding cultures")
            end
            age = age(age_count > 0);
            
            culture_table_array = {};
            culture_array = {};
            for iC = 1:length(rg.Cultures)
                rec_table_array = cell(1,length(age));
                culture_metadata = [rg.Cultures{iC}.Metadata];
                div = [culture_metadata.DIV];
                if length(div) >= length(age)
                    age_dev = abs(div - age');
                    age_count = sum(age_dev<=tolerance);
                    if sum(age_count == 1) == length(age) % Need to find a way to handle two recordings falling within the tolerance window
                        sel_rec = rg.Cultures{iC}(age_count == 1);
                        culture_array = [culture_array {sel_rec}];
                        for iR = 1:length(sel_rec)
                            iR_table = getFeatureTable(rg.Cultures{iC}(iR), unit_features, network_features);
                            iR_table.Properties.VariableNames = iR_table.Properties.VariableNames + "_" + string(div(iR));
                            rec_table_array{iR} = iR_table;
                        end
                        culture_table_array{iC} = [rec_table_array{:}];
                        
                    else
                        continue
                    end
                else
                   continue
               end
            end
            feature_table = vertcat(culture_table_array{:});
        end
        
        function [norm_train_data, norm_test_data] = normalizeByGroup(rg, feat_mat, object_group, grouping_var, train_idx, test_idx)
            arguments
                rg RecordingGroup
                feat_mat {isnumeric}
                object_group %Units/recordings/cultures corresponding to the rows of feat_mat
                grouping_var string = "PlatingDate" %Has to correspond to a MEArecording metadata field
                train_idx (1,:) {isnumeric} = ones(1,size(feat_mat,1)) %Default to unsupervised // should be binary and length of object_group
                test_idx (1,:) {isnumeric} = zeros(1,size(feat_mat,1))
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
                iBatch_train = (iG == g) & train_idx;
                batch_train_data = feat_mat(iBatch_train,:);
                batch_mean = mean(batch_train_data);
                batch_std = std(batch_train_data);
                feat_mat(iBatch_train,:) = (batch_train_data - batch_mean) ./ batch_std;
                if any(test_idx)
                    iBatch_test = iG == g & test_idx;
                    batch_test_data = feat_mat(iBatch_test,:);
                    feat_mat(iBatch_test,:) = (batch_test_data - batch_mean) ./ batch_std;
                end
            end
            if length(G) > 1 %Only normalize again if more than 1 group exists
                train_data = feat_mat(train_idx,:);
                train_mean = mean(train_data);
                train_std = std(train_data);
                norm_train_data = (train_data - train_mean) ./ train_std;
                test_data = feat_mat(test_idx,:);
                norm_test_data = (test_data - train_mean) ./ train_std;
            else
                norm_train_data = feat_mat(train_idx,:);
                norm_test_data = feat_mat(test_idx,:);
            end
        end
        
        function reduction = reduceDimensionality(rg, level, method, n_dims, grouping_var, unit_features, network_features, age, tolerance)
            arguments
               rg RecordingGroup
               level string = "Unit" %Unit, Recording or Culture level
               method string = "UMAP" %UMAP, PCA
               n_dims (1,1) {isnumeric} = 2 %Number of output dimensions
               grouping_var string = "PlatingDate" %Has to correspond to a MEArecording metadata field
               unit_features string = "all"
               network_features string = "all"
               age {isnumeric} = [7,14,21,28] %Only relevant for DimensionalityReduction on the culture level
               tolerance {isnumeric} = 1 %Gives tolerance for culture selection by age (e.g. age=7 tolerance=1 allows DIVs 6-8)
            end
            
            fprintf('Performing %s dimensionality reduction on %ss and normalizing by %s\n', method, level, grouping_var)
            
            switch level
                case "Unit"
                    object_group = [rg.Units];
                    ref_wf = [object_group.ReferenceWaveform];
                    ref_wf = ref_wf(sum(ref_wf,2)~=0,:);
                    unit_feature_array = [object_group.AverageFeatures];
                    if isempty(unit_feature_array)
                        fprintf('No unit features found. Continuing using only the reference waveform. \n')
                        input_mat = double(ref_wf');
                    else
                        feature_table = rg.Recordings.concatenateFeatures(unit_feature_array,"Activity");
                        aligned_wf = rg.alignWaveforms();
                        input_mat = double([aligned_wf' feature_table.Variables]);
%                         feature_table = rg.Recordings.concatenateFeatures(unit_feature_array,"all");
%                         input_mat = feature_table.Variables;
                    end
                    
                case "Recording"
                    object_group = rg.Recordings;
                    feature_table = rg.aggregateRecordingFeatureTables(unit_features, network_features);
                    input_mat = feature_table.Variables;
                    
                case "Culture" 
                    [feature_table, object_group] = aggregateCultureFeatureTables(rg, age, tolerance, unit_features, network_features);
                    input_mat = feature_table.Variables;
            end
            
            input_mat(isnan(input_mat)) = 0;%Handle rare case where NaN appears
            norm_data = normalizeByGroup(rg, input_mat, object_group, grouping_var); %Normalize
            norm_data(:,any(isnan(norm_data))) = [];%Remove peak NaNs
            norm_data = norm_data./max(abs(norm_data)); %Scale data
            
            switch method
                case "UMAP"
                    color_file = fullfile(rg.Recordings(1).getParentPath(),'umap','colorsByName.properties');
                    [reduction, umap, clusterIdentifiers, extras] = run_umap(norm_data,'n_components',n_dims,'n_neighbors',100,'min_dist',0.1,'cluster_detail','adaptive','spread',1,'sgd_tasks',20,...
                        'verbose','none','color_file',color_file);
                    
                case "PCA"
                    [coeff,reduction,latent] = pca(norm_data);
                    reduction = reduction(:,1:n_dims);
                    
            end
            
            rg.DimensionalityReduction.(level).(method).Reduction = reduction;
            rg.DimensionalityReduction.(level).(method).GroupingVariable = grouping_var;
            rg.DimensionalityReduction.(level).(method).UnitFeatures = unit_features;
            rg.DimensionalityReduction.(level).(method).NetworkFeatures = network_features;
            rg.DimensionalityReduction.(level).ObjectGroup = object_group;
            rg.plot_dimensionality_reduction(reduction);
            %             cluster_idx = num2cell(cluster_idx); %Prepare to use deal to assign cluster ids
            %             [rg.Units.ClusterID] = deal(cluster_idx{:});
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
               metadata_object = [rg.Units.MEArecording];
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
                [idx_array{m}, group_array{m}] = returnMetadataArray(rg, metadata_object, metadata_names(m));
            end
            group_vals = cellfun(@unique, idx_array, "un",0);
            group_combs = combvec(group_vals{:});
            group_idx = vertcat(idx_array{:});
            [~,cluster_idx] = ismember(group_idx', group_combs','rows');
            group_labels = arrayfun(@(x) group_array{x}(group_combs(x,:))',1:length(group_array),"un",0);
            group_labels_comb = join([group_labels{:}]);
            
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
        
        function aligned_wf = alignWaveforms(rg)
            object_group = [rg.Units];
            ref_wf = [object_group.ReferenceWaveform];
            ref_wf = ref_wf(sum(ref_wf,2)~=0,:);
            [~,i] = min(ref_wf,[],1);
            peak_idx = mean(i);
            max_offset = round(peak_idx/2);
            x = max_offset:size(ref_wf,1)+max_offset-1;
            xq = 1:size(ref_wf,1)+2*max_offset;
            
            interp_wf = interp1(x,ref_wf,xq,"linear",'extrap');
            rm_idx = find(abs(peak_idx - i) >= max_offset);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plot_true_clusters(rg, level, method, grouping_var)
            arguments
                rg RecordingGroup
                level string = "Unit" %Unit, Recording, Culture
                method string = "UMAP" %dim reduction method: "UMAP","PCA"
                grouping_var string = "CellLine"
            end
            reduction = rg.DimensionalityReduction.(level).(method).Reduction;
            metadata_object = rg.DimensionalityReduction.(level).ObjectGroup;
            [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, grouping_var);
            plot_cluster_outlines(rg,reduction, cluster_idx)
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
               cmap (:,3) {isnumeric} = othercolor('Set19',length(unique(cluster_idx)));
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
        
        function plot_dimensionality_reduction(rg,reduction,cluster_idx)
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
            figure('Color','w');
            if size(reduction,2) == 2
                nexttile
                scatter(reduction(:,1),reduction(:,2),sz,cluster_idx,'filled')
                c_map = othercolor('Set19',N_clust);
                colormap(c_map)
                cb = colorbar; cb.Title.String = "Cluster Index";
                cb.Ticks = linspace(cb.Limits(1) + (N_clust-1)/(2*N_clust), cb.Limits(2) - (N_clust-1)/(2*N_clust), N_clust); cb.TickLabels = 1:N_clust;
                if N_clust > 1
                    for i = 1:N_clust
                        nexttile
                        scatter(reduction(cluster_idx==i,1),reduction(cluster_idx==i,2),sz,'filled')
                    end
                end
            else
                nexttile
                scatter3(reduction(:,1),reduction(:,2),reduction(:,3),sz,'filled')
                c_map = othercolor('Set19',N_clust);
                colormap(c_map)
                cb = colorbar; cb.Title.String = "Cluster Index";
                cb.Ticks = linspace(cb.Limits(1) + (N_clust-1)/(2*N_clust), cb.Limits(2) - (N_clust-1)/(2*N_clust), N_clust); cb.TickLabels = 1:N_clust;
                if N_clust > 1
                    for i = 1:N_clust
                        nexttile
                        scatter3(reduction(cluster_idx==i,1),reduction(cluster_idx==i,2),reduction(cluster_idx==i,3),sz,'filled')
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
            aligned_wf = rg.alignWaveforms();
            N_clust = length(unique(cluster_idx));
            avg_wf = cell(1,N_clust);
            cluster_wfs = cell(1,N_clust);
            cluster_size = groupcounts(cluster_idx);
            time_vector = 1:size(aligned_wf,1)/20;
            
            for c = 1:N_clust
                cluster_wfs{c} = aligned_wf(:,cluster_idx == c);
                avg_wf{c} = median(cluster_wfs{c},2);
            end
            
            figure("Color","w");
            nexttile
            plot(time_vector,horzcat(avg_wf{:}),"LineWidth",2)
            legend(num2str(cluster_size))
            box off; xlabel("Time [ms]")
            for c = 1:N_clust
                nexttile
                p = plot(time_vector,cluster_wfs{c},"k","LineWidth",0.1,'Color',[0 0 0 0.1]);
                hold on
                plot(time_vector,horzcat(avg_wf{c}),"LineWidth",2)
                box off; xlabel("Time [ms]")
            end
        end
    end
end