classdef RecordingGroup < handle
    properties
        Recordings
        Units
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
            rg.Recordings = rg.filterRecordingArray(recording_array);
            rg.Units = [rg.Recordings.Units];
            fprintf('Initialized RecordingGroup with %i recordings and %i units\n',...
                length(rg.Recordings),length(rg.Units))
        end
        
        function parseParameters(obj, parameters)
            obj.Parameters = obj.returnDefaultParams();
            parameter_fields = fieldnames(parameters);
            if isempty(parameter_fields)
                warning("No parameters provided, using default ones")
            else
                for pf = parameter_fields
                    if ismember(pf,fieldnames(obj.Parameters))
                        subfields = fieldnames(parameters.(pf));
                        for sf = subfields
                            if ismember(sf,fieldnames(obj.Parameters.(pf)))
                                obj.Parameters.(pf).(sf) = parameters.(pf).(sf);
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
        
        function array = filterRecordingArray(obj, recording_array)
            % Input as cell arrays
            inclusion = obj.Parameters.Selection.Inclusion;
            for i = 1:length(inclusion)
                if isnumeric(inclusion{i}{2})
                    value_vector = [recording_array.(inclusion{i}{1})];
                    idx = ismember(value_vector,[inclusion{i}{2:end}]);
                else
                    value_vector = {recording_array.(inclusion{i}{1})};
                    idx = ismember(value_vector,inclusion{i}(2:end));
                end
                recording_array = recording_array(idx);
            end
            
            exclusion = obj.Parameters.Selection.Exclusion;
            for i = 1:length(exclusion)
                if isnumeric(exclusion{i}{2})
                    value_vector = [recording_array.(exclusion{i}{1})];
                    idx = ismember(value_vector,[exclusion{i}{2:end}]);
                else
                    value_vector = {recording_array.(exclusion{i}{1})};
                    idx = ismember(value_vector,exclusion{i}(2:end));
                end
                recording_array = recording_array(~idx);
            end
            array = recording_array;
        end
    end
    
    methods (Static)
       function defaultParams = returnDefaultParams()
           %Input to the filterRecordingArray function
           defaultParams.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value
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
               normalization string = "PlatingDate" %Has to correspond to MEArecording metadata field
            end
            
            fprintf('Performing %s dimensionality reduction on %ss and normalizing by %s\n', method, level, normalization)
            
            switch level
                case "unit"
                    unit_array = [rg.Units];
                    ref_wf = [unit_array.ReferenceWaveform];
                    unit_feature_array = [unit_array.Features];
                    input_mat = rg.Recordings.concatenateFeatures(unit_feature_array,["Activity"]);
                    norm_mat = normalize(input_mat.Variables,1);
                    input_mat = double([ref_wf' norm_mat]);
                case "recording"
                    
                case "culture" %Create function to group recordings of same culture
                    
            end
            switch method
                case "umap"
                    [reduction, umap, clusterIdentifiers, extras]=run_umap(input_mat,'n_components',2,'n_neighbors',100,'min_dist',0.1,'cluster_detail','adaptive','spread',1,'sgd_tasks',20,...
                        'verbose','none','color_file','C:\Users\Philipp\Downloads\umapAndEppFileExchange (4.1)\colorsByName.properties');
                    cluster_idx = kmeans(reduction,2);
                case "pca"
                    [coeff,reduction,latent] = pca(input_mat);
                    cluster_idx = kmeans(reduction,2);
            end
        end
    end
end