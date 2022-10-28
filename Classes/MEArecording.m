classdef MEArecording < handle
    properties
        Metadata
        RecordingInfo
        Parameters
        Units
        NetworkFeatures
        UnitFeatures
        Connectivity
        GraphFeatures
    end
    
    methods (Hidden)
        function mearec = MEArecording(metadata, parameters, analyses)
            arguments
                metadata (1,1) struct %Metadata to be parsed // InputPath is required
                parameters (1,1) struct = struct()
                analyses (1,1) struct = struct() %Implement way to toggle analyses on/off
            end
            
            if nargin > 0
                mearec.parseMetadata(metadata);
                mearec.parseRecordingInfo();
                mearec.parseParameters(parameters);
                mearec.performAnalyses(analyses);
            end
        end
        
        function parseMetadata(obj, metadata)
            if isempty(fieldnames(metadata)) %Check if any metadata is provided
                error("No metadata provided, cannot continue");
            elseif ~isfield(metadata,"InputPath")
                error("Metadata does not provide InputPath, cannot continue");
            else
                obj.Metadata = metadata;
                disp("Imported metadata")
            end
        end
        
        function parseRecordingInfo(obj)
            obj.RecordingInfo.ElectrodeCoordinates = obj.getElectrodeCoordinates();
            obj.RecordingInfo.SamplingRate = obj.getSamplingRate();
            obj.RecordingInfo.Duration = obj.getRecordingDuration();
            disp("Loaded recording information")
        end
        
        function xy_coordinates = getElectrodeCoordinates(obj)
            xy_coordinates = readNPY(fullfile(obj.Metadata.InputPath, "channel_positions.npy"));
        end
        
        function sampling_rate = getSamplingRate(obj)
            params_path = fullfile(obj.Metadata.InputPath, "params.py");
            fileID = fopen(params_path,'r');
            params_txt = fscanf(fileID,'%c');
            params_parts = strsplit(params_txt,'\n');
            sr_idx = startsWith(params_parts,'sample');
            sr = strsplit(params_parts{sr_idx},' = ');
            sampling_rate = str2double(sr{2});
        end
        
        function duration = getRecordingDuration(obj)
            spike_times = double(readNPY(fullfile(obj.Metadata.InputPath, "spike_times.npy")));
            duration = ceil(max(spike_times)/obj.RecordingInfo.SamplingRate);
        end
        
        function parseParameters(obj, parameters)
            obj.Parameters = obj.returnDefaultParams();
            parameter_fields = fieldnames(parameters);
            if isempty(parameter_fields)
                warning("No parameters provided, using default ones")
            else
                for pf = parameter_fields
                    subfields = fieldnames(parameters.(pf));
                    for sf = subfields
                        obj.Parameters.(pf).(sf) = parameters.(pf).(sf);
                    end
                end
                disp("Imported custom parameters")
            end
        end
        
        function performAnalyses(obj, analyses)
            obj.generateUnits()
        end
        
        function importUnits(obj)
            [max_amplitude, norm_wf_matrix] = obj.generateWaveformMatrix();
            
        end
        
        function [max_amplitude, norm_wf_matrix] = generateWaveformMatrix(obj)
            template_matrix = readNPY(fullfile(obj.Metadata.InputPath, "templates.npy"));
            [max_amplitude, max_idx] = min(template_matrix,[],[2,3],'linear');
            max_amplitude = max_amplitude * obj.Parameters.QC.lsb;
            [~,~,reference_electrode] = ind2sub(size(template_matrix),max_idx);
            peak_wf = arrayfun(@(x) template_matrix(x,:,reference_electrode(x)),1:length(reference_electrode),'un',0);
            peak_wf_matrix = cat(1,peak_wf{:})';
            norm_wf_matrix = peak_wf_matrix./max(abs(peak_wf_matrix));
        end
    end
    
    methods (Static)
        function defaultParams = returnDefaultParams()
            % Parameters for the unit quality control (QC)
            defaultParams.QC.lsb = 6.2; %Factor by which to multiply voltage data
            defaultParams.QC.amp = NaN; %[Minimum Maximum] amplitude
            defaultParams.QC.rate = [0.1 10]; %[Minimum Maximum] firing rate
            defaultParams.QC.rv = 0.02; %Refractory period violation ratio
            
            % Parameters for the burst detection
            defaultParams.Bursts.N = 0.0015; %Number of spikes in bursts
            defaultParams.Bursts.ISI_N = 1.5; %Burst length
            defaultParams.Bursts.merge_t = 2; %IBI to be merged
            defaultParams.Bursts.binning = 0.1;
            
            % Parameters for the burst detection
            defaultParams.Regularity.binning = 0.1; %Binning to compute regularity features
        end
    end
end