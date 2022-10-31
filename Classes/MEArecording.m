classdef MEArecording < handle
    properties
        Metadata
        RecordingInfo
        Parameters
        Spikes
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
            [obj.Spikes.Times, obj.Spikes.Units] = obj.getSpikeTimes();
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
        
        function [spike_times, spike_units] = getSpikeTimes(obj)
            spike_times = double(readNPY(fullfile(obj.Metadata.InputPath, "spike_times.npy"))) / obj.RecordingInfo.SamplingRate;
            spike_units = double(readNPY(fullfile(obj.Metadata.InputPath, "spike_templates.npy")))+1;
        end
        
        function template_matrix = getTemplateMatrix(obj)
            template_matrix = readNPY(fullfile(obj.Metadata.InputPath, "templates.npy"));
        end
        
        function duration = getRecordingDuration(obj)
            duration = ceil(max(obj.Spikes.Times));
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
        
        function performAnalyses(obj, analyses)
            if analyses.SingleCell
                obj.Units = obj.generateUnits();
                obj.updateSpikeTimes();
            end
            
            if analyses.Regularity
                
            end
        end
        
        function unit_array = generateUnits(obj)
            [max_amplitudes, norm_wf_matrix] = obj.generateWaveformMatrix();
            [firing_rates, unit_spike_times] = obj.calculateFiringRates();
            good_units = performUnitQC(obj, max_amplitudes, firing_rates, unit_spike_times, norm_wf_matrix);
            if sum(good_units) > 0
                good_amplitudes = max_amplitudes(good_units);
                good_unit_spike_times = unit_spike_times(good_units);
                good_wf_matrix = norm_wf_matrix(:,good_units);
                waveform_features = obj.inferWaveformFeatures(good_amplitudes, good_wf_matrix);
                unit_array = Unit();
                for u = 1:length(waveform_features)
                    unit_array(u) = Unit(obj, good_wf_matrix(:,u), good_unit_spike_times{u}, waveform_features{u});
                end
            else
                warning("No good units found")
            end
        end
        
        function [max_amplitudes, norm_wf_matrix] = generateWaveformMatrix(obj)
            template_matrix = obj.getTemplateMatrix();
            [max_amplitudes, max_idx] = min(template_matrix,[],[2,3],'linear');
            max_amplitudes = abs(max_amplitudes * obj.Parameters.QC.LSB);
            [~,~,reference_electrode] = ind2sub(size(template_matrix),max_idx);
            peak_wf = arrayfun(@(x) template_matrix(x,:,reference_electrode(x)),1:length(reference_electrode),'un',0);
            peak_wf_matrix = cat(1,peak_wf{:})';
            norm_wf_matrix = peak_wf_matrix./max(abs(peak_wf_matrix));
        end
        
        function [firing_rates, unit_spike_times] = calculateFiringRates(obj)
            unit_spike_times = splitapply(@(x) {x}, obj.Spikes.Times, obj.Spikes.Units);
            firing_rates = cellfun(@length, unit_spike_times)/obj.RecordingInfo.Duration;
        end
        
        function good_units = performUnitQC(obj, max_amplitudes, firing_rates, unit_spike_times, norm_wf_matrix)
            bad_amplitude = checkAmplitude(obj, max_amplitudes);
            bad_firing_rate = checkFiringRate(obj,firing_rates);
            bad_rpv = checkRefractoryPeriodViolations(obj, unit_spike_times);
            [is_axonal, is_noisy] = checkWaveform(obj, norm_wf_matrix);
            good_units = ~(bad_amplitude | bad_firing_rate | bad_rpv | is_axonal | is_noisy);
            fprintf('Found %i/%i good units\n',sum(good_units),length(good_units))
        end
        
        function bad_amplitude = checkAmplitude(obj, max_amplitudes)
            if ~isempty(obj.Parameters.QC.Amplitude) && all(~isnan(obj.Parameters.QC.Amplitude))
                bad_amplitude = max_amplitudes < obj.Parameters.QC.Amplitude(1) | max_amplitudes > obj.Parameters.QC.Amplitude(2);
            else
                bad_amplitude = zeros(size(max_amplitudes));
            end
            fprintf('Found %i units with bad amplitudes\n',sum(bad_amplitude))
        end
        
        function bad_firing_rate = checkFiringRate(obj, firing_rates)
            if ~isempty(obj.Parameters.QC.FiringRate) && all(~isnan(obj.Parameters.QC.FiringRate))
                bad_firing_rate = firing_rates < obj.Parameters.QC.FiringRate(1)| firing_rates > obj.Parameters.QC.FiringRate(2);
            else
                bad_firing_rate = zeros(size(firing_rates));
            end
            fprintf('Found %i units with bad firing rates\n',sum(bad_firing_rate))
        end
        
        function bad_rpv = checkRefractoryPeriodViolations(obj, unit_spike_times)
            if ~isempty(obj.Parameters.QC.RPV) && ~isnan(obj.Parameters.QC.RPV)
                rpv_ratio = cellfun(@(x) sum(diff(x)<0.002)/(length(x)-1),unit_spike_times); %2 ms RP
                bad_rpv = rpv_ratio > obj.Parameters.QC.RPV;
            else
                bad_rpv = zeros(size(unit_spike_times));
            end
            fprintf('Found %i units with bad RPV\n',sum(bad_rpv))
        end
        
        function [is_axonal, is_noisy] = checkWaveform(obj, norm_wf_matrix)
            if ~isnan(obj.Parameters.QC.Axon)
                amplitude_ratio = -(max(norm_wf_matrix)./min(norm_wf_matrix))';
                %Fit two gaussians and try to find cutoff
                is_axonal = amplitude_ratio > obj.Parameters.QC.Axon;
            else
                is_axonal = zeros(size(norm_wf_matrix,2),1);
            end
            fprintf('Identified %i units as axonal\n',sum(is_axonal))
            
            if ~isnan(obj.Parameters.QC.Noise)
                [~,peak_amp_idx] = min(norm_wf_matrix);
                cutout_idx = peak_amp_idx + obj.Parameters.QC.NoiseCutout(1) * obj.RecordingInfo.SamplingRate/1000:...
                    peak_amp_idx + obj.Parameters.QC.NoiseCutout(2) * obj.RecordingInfo.SamplingRate/1000;
                noise_indicator = sum(abs(diff(diff(norm_wf_matrix(cutout_idx,:))>0)))'; %Check for non-monotomies in waveform shapes as noise indicators
                
                is_noisy = noise_indicator > obj.Parameters.QC.Noise;
            else
                is_noisy = zeros(size(norm_wf_matrix,2),1);
            end
            fprintf('Identified %i units as noise\n',sum(is_noisy))
        end
        
        function waveform_features = inferWaveformFeatures(obj,max_amplitudes, norm_wf_matrix)
            interpolation_factor = 10;
            ms_conversion = (obj.RecordingInfo.SamplingRate / 1000) * interpolation_factor;
            x = 1:size(norm_wf_matrix,1);
            xq = 1:(1/interpolation_factor):size(norm_wf_matrix,1);
            interp_wf_matrix = double(interp1(x,norm_wf_matrix,xq,'pchip'));
            zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %Function to detect zero crossings
            tx = zci(interp_wf_matrix); %Apply to interpolated waveform matrix
            [sample_idx,electrode_idx] = ind2sub(size(interp_wf_matrix),tx);
            unit_zero_crossings = splitapply(@(x) {x},sample_idx,electrode_idx);
            [unit_trough_value,unit_trough_idx] = min(interp_wf_matrix);
            peak_1_cutout = interp_wf_matrix(unit_trough_idx - ms_conversion:unit_trough_idx,:); %Limit detection range for peaks
            peak_2_cutout = interp_wf_matrix(unit_trough_idx:unit_trough_idx + ms_conversion,:);
            [peak_1_value, ~] = max(peak_1_cutout);
            [peak_2_value, peak_2_idx] = max(peak_2_cutout);
            peak_2_idx = peak_2_idx + unit_trough_idx - 1;
            asymmetry = (peak_2_value - peak_1_value) ./ (peak_2_value + peak_1_value);
            t2pdelay = (peak_2_idx - unit_trough_idx) ./ (obj.RecordingInfo.SamplingRate / 1000); %in [ms]
            t2pratio = abs(unit_trough_value ./ peak_2_value);
            waveform_features = cell(1,length(max_amplitudes));
            for u = 1:size(interp_wf_matrix,2)
                zc = [1 unit_zero_crossings{u}' length(xq)];
                [~,zc_pre_trough] = max(zc(zc<unit_trough_idx(u))); %Find zero crossing towards the trough
                unit_features.AUC_peak_1 = trapz(interp_wf_matrix(zc(zc_pre_trough - 1):zc(zc_pre_trough)));
                unit_features.AUC_trough = abs(trapz(interp_wf_matrix(zc(zc_pre_trough):zc(zc_pre_trough + 1))));
                unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):zc(zc_pre_trough + 2)));
                unit_features.Rise = slewrate(interp_wf_matrix(unit_trough_idx(u):peak_2_idx(u),u));
                unit_features.Decay = slewrate(interp_wf_matrix(peak_2_idx(u):end,u));
                unit_features.Asymmetry = asymmetry(u);
                unit_features.T2Pdelay = t2pdelay(u);
                unit_features.T2Pratio = t2pratio(u);
                waveform_features{u} = unit_features;
            end
        end
        
        function updateSpikeTimes(obj)
           % Update spike times and units after the quality control 
           spike_times = vertcat(obj.Units.SpikeTimes);
           spike_units = arrayfun(@(x) ones(1,length(obj.Units(x).SpikeTimes))*x,1:length(obj.Units),'un',0);
           spike_units = horzcat(spike_units{:});
           [spike_times_sorted, sort_idx] = sort(spike_times,'ascend');
           spike_units_sorted = spike_units(sort_idx);
           obj.Spikes.Times = spike_times_sorted;
           obj.Spikes.Units = spike_units_sorted';
        end
        
        function getRegularity(obj)
            
        end
    end
    
    methods (Static)
        function defaultParams = returnDefaultParams()
            % Parameters for the unit quality control (QC)
            defaultParams.QC.LSB = 6.2; %Factor by which to multiply voltage data
            defaultParams.QC.Amplitude = [20 300]; %[Minimum Maximum] amplitude
            defaultParams.QC.FiringRate = [0.1 10]; %[Minimum Maximum] firing rate
            defaultParams.QC.RPV = 0.02; %Refractory period violation ratio
            defaultParams.QC.Axon = 0.8; %Ratio (positive amplitude)/ (maximum amplitude) to be considered axonal
            defaultParams.QC.Noise = 14; % # of sign changes for the waveform derivative ("up-and-down" waveform)
            defaultParams.QC.NoiseCutout = [-1 1]; % [ms] before and after peak to be considered for noise detection
            
            % Parameters for the burst detection
            defaultParams.Bursts.N = 0.0015; %Number of spikes in bursts
            defaultParams.Bursts.ISI_N = 1.5; %Burst length
            defaultParams.Bursts.merge_t = 2; %IBI to be merged
            defaultParams.Bursts.binning = 0.1;
            
            % Parameters for the burst detection
            defaultParams.Regularity.binning = 0.1; %Binning to compute regularity features
        end
    end
    
    methods
       % Plots
       function NetworkScatterPlot(obj)
          figure('Color','w');
          plot(obj.Spikes.Times,obj.Spikes.Units,'k.')
          xlabel('Time [s]')
          ylabel('Unit ID')
          axis tight
          box off
       end
    end
end