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
        
        function generateUnits(obj)
            [max_amplitude, norm_wf_matrix] = obj.generateWaveformMatrix();
            [firing_rates, unit_spike_times] = obj.calculateFiringRates();
            performUnitQC(obj, max_amplitude, firing_rates, unit_spike_times, norm_wf_matrix)
        end
        
        function [max_amplitudes, norm_wf_matrix] = generateWaveformMatrix(obj)
            template_matrix = readNPY(fullfile(obj.Metadata.InputPath, "templates.npy"));
            [max_amplitudes, max_idx] = min(template_matrix,[],[2,3],'linear');
            max_amplitudes = abs(max_amplitudes * obj.Parameters.QC.LSB);
            [~,~,reference_electrode] = ind2sub(size(template_matrix),max_idx);
            peak_wf = arrayfun(@(x) template_matrix(x,:,reference_electrode(x)),1:length(reference_electrode),'un',0);
            peak_wf_matrix = cat(1,peak_wf{:})';
            norm_wf_matrix = peak_wf_matrix./max(abs(peak_wf_matrix));
        end
        
        function [firing_rates, unit_spike_times] = calculateFiringRates(obj)
            spike_times = double(readNPY(fullfile(obj.Metadata.InputPath, "spike_times.npy")))/obj.RecordingInfo.SamplingRate;
            spike_unit = double(readNPY(fullfile(obj.Metadata.InputPath, "spike_templates.npy")))+1;
            unit_spike_times = splitapply(@(x) {x}, spike_times, spike_unit);
            firing_rates = cellfun(@length, unit_spike_times)/obj.RecordingInfo.Duration;
        end
        
        function performUnitQC(obj, max_amplitude, firing_rates, unit_spike_times, norm_wf_matrix)
            bad_amplitude = checkAmplitude(obj, max_amplitude);
            bad_firing_rate = checkFiringRate(obj,firing_rates);
            bad_rpv = checkRefractoryPeriodViolations(obj, unit_spike_times);
            [is_axonal, is_noisy] = checkWaveform(obj, norm_wf_matrix);
        end
        
        function bad_amplitude = checkAmplitude(obj, max_amplitudes)
            if ~isempty(obj.Parameters.QC.Amplitude) && all(~isnan(obj.Parameters.QC.Amplitude))
                bad_amplitude = max_amplitudes > obj.Parameters.QC.Amplitude(1) & max_amplitudes < obj.Parameters.QC.Amplitude(2);
            else
                bad_amplitude = ones(size(max_amplitudes));
            end
        end
        
        function bad_firing_rate = checkFiringRate(obj, firing_rates)
            if ~isempty(obj.Parameters.QC.FiringRate) && all(~isnan(obj.Parameters.QC.FiringRate))
                bad_firing_rate = firing_rates > obj.Parameters.QC.FiringRate(1) & firing_rates < obj.Parameters.QC.FiringRate(2);
            else
                bad_firing_rate = ones(size(firing_rates));
            end
        end
        
        function bad_rpv = checkRefractoryPeriodViolations(obj, unit_spike_times)
            if ~isempty(obj.Parameters.QC.RPV) && ~isnan(obj.Parameters.QC.RPV)
                rpv_ratio = cellfun(@(x) sum(diff(x)<0.002)/(length(x)-1),unit_spike_times); %2 ms RP
                bad_rpv = rpv_ratio > obj.Parameters.QC.RPV;
            else
                bad_rpv = ones(size(unit_spike_times));
            end
        end
        
        function [is_axonal, is_noisy] = checkWaveform(obj, norm_wf_matrix)
            if ~isnan(obj.defaultParams.QC.Axon)
                amplitude_ratio = -(max(norm_wf_matrix)./min(norm_wf_matrix));
                %Fit two gaussians and try to find cutoff
            else
                
            end
            
            if ~isnan(defaultParams.QC.Noise)
                noise_indicator = sum(abs(diff(diff(norm_wf_matrix(:,15:50)')>0))); %Check for non-monotomies in waveform shapes as noise indicators
                noise_idx = noise_indicator>noise_cutoff;
            else
                
            end
        end
    end
    methods (Static)
        function defaultParams = returnDefaultParams()
            % Parameters for the unit quality control (QC)
            defaultParams.QC.LSB = 6.2; %Factor by which to multiply voltage data
            defaultParams.QC.Amplitude = NaN; %[Minimum Maximum] amplitude
            defaultParams.QC.FiringRate = [0.1 10]; %[Minimum Maximum] firing rate
            defaultParams.QC.RPV = 0.02; %Refractory period violation ratio
            defaultParams.QC.Axon = 0.8; %Ratio (positive amplitude)/ (maximum amplitude) to be considered axonal
            defaultParams.QC.Noise = 14; % # of sign changes for the waveform derivative ("up-and-down" waveform)
            
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