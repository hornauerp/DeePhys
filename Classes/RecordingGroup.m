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
            rg.Cultures = rg.groupCultures();
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
                if isnumeric(inclusion{i}{2})
                    value_vector = [metadata_array.(inclusion{i}{1})];
                    idx = ismember(value_vector,[inclusion{i}{2:end}]);
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
                if isnumeric(exclusion{i}{2})
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
            chip_ids = unique([metadata.ChipID]);
            for id = chip_ids
                inclusion = {{'ChipID',id}};
                chip_idx = rg.filterRecordingArray(metadata, inclusion);
                chip_array = rg.Recordings(chip_idx);
                chip_metadata = [chip_array.Metadata];
                chip_plating_dates = unique([chip_metadata.PlatingDate]);
                for pd = chip_plating_dates
                   inclusion = {{'PlatingDate',pd}};
                   plating_idx = rg.filterRecordingArray(metadata, inclusion);
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
        function feature_table = aggregateFeatureTables(rg, unit_features, network_features)
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
        
        function [reduction, cluster_idx] = reduceDimensionality(rg, level, method, normalization)
            arguments
               rg RecordingGroup
               level string = "unit" %Unit, recording or culture level
               method string = "umap" %UMAP, PCA
               normalization string = "PlatingDate" %Has to correspond to a MEArecording metadata field
            end
            
            fprintf('Performing %s dimensionality reduction on %ss and normalizing by %s\n', method, level, normalization)
            
            switch level
                case {"unit","Unit","u","U"}
                    unit_array = [rg.Units];
                    ref_wf = [unit_array.ReferenceWaveform];
                    unit_feature_array = [unit_array.Features];
                    if isempty(unit_feature_array)
                        fprintf('No unit features found. Continuing using only the reference waveform. \n')
                        input_mat = double(ref_wf');
                    else
                        input_mat = rg.Recordings.concatenateFeatures(unit_feature_array,"Activity");
                        norm_mat = normalize(input_mat.Variables,1);
                        input_mat = double([ref_wf' norm_mat]);
                    end
                    
                case {"recording","Recording","r","R"}
                    
                case {"culture","Culture","c","C"} %Create function to group recordings of same culture
                    
            end
            
            switch method
                case {"umap","UMAP","Umap"}
                    color_file = fullfile(rg.Recordings(1).getParentPath(),'umap','colorsByName.properties');
                    [reduction, umap, clusterIdentifiers, extras] = run_umap(input_mat,'n_components',2,'n_neighbors',100,'min_dist',0.1,'cluster_detail','adaptive','spread',1,'sgd_tasks',20,...
                        'verbose','none','color_file',color_file);
                    cluster_idx = clusterIdentifiers;%kmeans(reduction,2);
                case {"pca","PCA","Pca"}
                    [coeff,reduction,latent] = pca(input_mat);
                    cluster_idx = kmeans(reduction,2);
            end
            rg.plot_cluster_scatter(reduction,cluster_idx);
            cluster_idx = num2cell(cluster_idx); %Prepare to use deal to assign cluster ids
            
            switch level
                case {"unit","Unit","u","U"}
                    [rg.Units.ClusterID] = deal(cluster_idx{:});
                case {"recording","Recording","r","R"}
                    [rg.Recordings.ClusterID] = deal(cluster_idx{:});
                case {"culture","Culture","c","C"} %Create function to group recordings of same culture
                    
            end
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_cluster_scatter(rg,reduction,cluster_idx)
            figure('Color','w');
            if size(reduction,2) == 2
                scatter(reduction(:,1),reduction(:,2),5,cluster_idx,'filled')
            else
                scatter3(reduction(:,1),reduction(:,2),reduction(:,3),5,cluster_idx,'filled')
            end
            N_clust = double(max(cluster_idx));
            c_map = othercolor('Set19',N_clust);
            colormap(c_map)
            cb = colorbar; cb.Title.String = "Cluster Index";
            cb.Ticks = linspace(cb.Limits(1) + (N_clust-1)/(2*N_clust), cb.Limits(2) - (N_clust-1)/(2*N_clust), N_clust); cb.TickLabels = 1:N_clust;
        end
    end
end