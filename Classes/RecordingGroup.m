classdef RecordingGroup < handle
    properties
        Recordings
        Units
        Cultures
        Parameters %Inclusion // exclusion criteria
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
            keep_recordings = rg.filterRecordingArray(recording_array);
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
                        subfields = fieldnames(parameters.(pf));
                        for sf = subfields
                            if ismember(sf,fieldnames(rg.Parameters.(pf)))
                                rg.Parameters.(pf).(sf) = parameters.(pf).(sf);
                            else
                                warning("Unrecognized parameter field: %s.%s",pf,sf)
                            end
                        end
                    else
                        warning("Unrecognized parameter field: %s",pf)
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
                chip_plating_dates = unique({chip_metadata.PlatingDate});
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
                case 'MEArecording'
                    metadata = [object_group.Metadata];
                case 'RecordingGroup'
                    metadata = cellfun(@(x) string(x(1).Metadata.(grouping_var)),object_group);
            end
            
            metadata = [metadata.(grouping_var)];
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
               age {isnumeric} = [7,14,21,28] %Only relevant for clustering on the culture level
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
                        input_mat = double([ref_wf' feature_table.Variables]);
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
            norm_data = norm_data./max(abs(norm_data)); %Scale data
            
            switch method
                case "UMAP"
                    color_file = fullfile(rg.Recordings(1).getParentPath(),'umap','colorsByName.properties');
                    [reduction, umap, clusterIdentifiers, extras] = run_umap(norm_data,'n_components',n_dims,'n_neighbors',100,'min_dist',0.1,'cluster_detail','adaptive','spread',1,'sgd_tasks',20,...
                        'verbose','none','color_file',color_file);
                    rg.Clustering.(level).Reduction.(method) = reduction;
                    
                case "PCA"
                    [coeff,reduction,latent] = pca(norm_data);
                    reduction = reduction(:,1:n_dims);
                    
            end
            rg.plot_cluster_scatter(reduction);
            %             cluster_idx = num2cell(cluster_idx); %Prepare to use deal to assign cluster ids
            %             [rg.Units.ClusterID] = deal(cluster_idx{:});
            
            rg.Clustering.(level).(method).Reduction = reduction;
            rg.Clustering.(level).(method).GroupingVariable = grouping_var;
            rg.Clustering.(level).(method).UnitFeatures = unit_features;
            rg.Clustering.(level).(method).NetworkFeatures = network_features;
            rg.Clustering.(level).ObjectGroup = object_group;
            
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
           end
           
           metadata_struct = [metadata_object.Metadata];
           assert(isfield(metadata_struct,metadata_name),"Not a valid metadata field")
           metadata = [metadata_struct.(metadata_name)];
           [iMetadata, metadata_groups] = findgroups(metadata);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plot_cluster_scatter(rg,reduction,cluster_idx)
            arguments
                rg RecordingGroup
                reduction {isnumeric} = rg.Clustering.reduction
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
    end
end