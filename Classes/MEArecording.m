classdef MEArecording < handle
    properties
        Metadata
        RecordingInfo
        Parameters
        Spikes
        Bursts
        Units
        NetworkFeatures
        UnitFeatures
        Connectivity
        GraphFeatures
    end
    
    methods (Hidden)
        function mearec = MEArecording(metadata, parameters)
            arguments
                metadata (1,1) struct %Metadata to be parsed // InputPath is required
                parameters (1,1) struct = struct()
            end
            
            if nargin > 0
                addpath(genpath(mearec.getParentPath()));
                mearec.parseMetadata(metadata);
                mearec.parseRecordingInfo();
                mearec.parseParameters(parameters);
                mearec.performAnalyses();
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
        
        function performAnalyses(obj)
            if obj.Parameters.Analyses.SingleCell
                obj.generateUnits();
                obj.aggregateSingleCellFeatures();
            end
            
            if obj.Parameters.Analyses.Regularity
                obj.getRegularity();
            end
            
            if obj.Parameters.Analyses.Bursts
                obj.getBurstStatistics();
            end
            
            if obj.Parameters.Analyses.Connectivity.CCG
                obj.inferConnectivityCCG();
            end
            
            if obj.Parameters.Analyses.Connectivity.DDC
                obj.inferConnectivityDDC();
            end
        end
       
        function [max_amplitudes, reference_electrode, norm_wf_matrix] = generateWaveformMatrix(obj)
            template_matrix = obj.getTemplateMatrix();
            [max_amplitudes, max_idx] = min(template_matrix,[],[2,3],'linear');
            max_amplitudes = abs(max_amplitudes * obj.Parameters.QC.LSB);
            [~,~,reference_electrode] = ind2sub(size(template_matrix),max_idx);
            peak_wf = arrayfun(@(x) template_matrix(x,:,reference_electrode(x)),1:length(reference_electrode),'un',0);
            peak_wf_matrix = cat(1,peak_wf{:})';
            norm_wf_matrix = peak_wf_matrix./max(abs(peak_wf_matrix));
        end
        
        function [firing_rates, unit_spike_times] = calculateFiringRates(obj)
            [spike_times, spike_units] = obj.getSpikeTimes();
            unit_spike_times = splitapply(@(x) {x}, spike_times, spike_units);
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
                unit_features = struct();
                zc = [1 unit_zero_crossings{u}' length(xq)];
                [~,zc_pre_trough] = max(zc(zc<unit_trough_idx(u))); %Find zero crossing towards the trough
                unit_features.AUC_peak_1 = trapz(interp_wf_matrix(zc(zc_pre_trough - 1):zc(zc_pre_trough)));
                unit_features.AUC_trough = abs(trapz(interp_wf_matrix(zc(zc_pre_trough):zc(zc_pre_trough + 1))));
                unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):zc(zc_pre_trough + 2)));
                %We pad the signal with min/max values to ensure reliable slewrate
                %calculations
                rise_cutout = interp_wf_matrix(unit_trough_idx(u):peak_2_idx(u),u);
                padded_rise = [ones(100,1)*rise_cutout(1); rise_cutout; ones(100,1)*rise_cutout(end)];
                unit_features.Rise = slewrate(padded_rise);
                decay_cutout = interp_wf_matrix(peak_2_idx(u):end,u);
                decay_cutout = decay_cutout(1:find(decay_cutout==min(decay_cutout)));
                padded_decay = [ones(100,1)*decay_cutout(1); decay_cutout; ones(100,1)*decay_cutout(end)];
                unit_features.Decay = slewrate(padded_decay,10);
                unit_features.Asymmetry = asymmetry(u);
                unit_features.T2Pdelay = t2pdelay(u);
                unit_features.T2Pratio = t2pratio(u);
                feature_names = fieldnames(unit_features);
                for f = 1:length(feature_names) %If slewrates could not be calculated // should happen very rarely
                   if isempty(unit_features.(feature_names{f}))
                       unit_features.(feature_names{f}) = 0;
                   end
                end
                waveform_features{u} = unit_features;
            end
        end
        
        function aggregateSingleCellFeatures(obj)
            feature_array = [obj.Units.Features];
            activity_array = [feature_array.ActivityFeatures];
            waveform_array = [feature_array.WaveformFeatures];
            activity_feature_names = fieldnames(obj.Units(1).Features.ActivityFeatures);
            waveform_feature_names = fieldnames(obj.Units(1).Features.WaveformFeatures);
            if ~isempty(obj.Parameters.Outlier.Method)
                for af = 1:length(activity_feature_names)
                    obj.UnitFeatures.Features.ActivityFeatures.(activity_feature_names{af}) = mean(rmoutliers([activity_array.(activity_feature_names{af})],...
                        obj.Parameters.Outlier.Method,'ThresholdFactor',obj.Parameters.Outlier.ThresholdFactor));
                end
                for wf = 1:length(waveform_feature_names)
                    obj.UnitFeatures.Features.WaveformFeatures.(waveform_feature_names{wf}) = mean(rmoutliers([waveform_array.(waveform_feature_names{wf})],...
                        obj.Parameters.Outlier.Method,'ThresholdFactor',obj.Parameters.Outlier.ThresholdFactor));
                end
            else
                for af = 1:length(activity_feature_names)
                    obj.UnitFeatures.Features.ActivityFeatures.(activity_feature_names{af}) = mean([activity_array.(activity_feature_names{af})]);
                end
                for wf = 1:length(waveform_feature_names)
                    obj.UnitFeatures.Features.WaveformFeatures.(waveform_feature_names{wf}) = mean([waveform_array.(waveform_feature_names{wf})]);
                end
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
        
        function Burst = detectBurstsISIN(obj, N, ISI_N)
            %Adapted from Bakkum et al., 2014
            tic;
            dT = zeros(N,length(obj.Spikes.Times))+inf;
            for j = 0:N-1
                dT(j+1,N:length(obj.Spikes.Times)-(N-1)) = obj.Spikes.Times( (N:end-(N-1))+j ) - ...
                    obj.Spikes.Times( (1:end-(N-1)*2)+j );
            end
            Criteria = zeros(size(obj.Spikes.Times)); % Initialize to zero
            Criteria( min(dT)<=ISI_N ) = 1; % Spike passes condition if it is
            % included in a set of N spikes
            % with ISI_N <= threshold.
            % %% Assign burst numbers to each spike
            SpikeBurstNumber = zeros(size(obj.Spikes.Times)) - 1; % Initialize to '-1'
            INBURST = 0; % In a burst (1) or not (0)
            NUM_ = 0; % Burst Number iterator
            NUMBER = -1; % Burst Number assigned
            BL = 0; % Burst Length
            for i = N:length(obj.Spikes.Times)
                if INBURST == 0 % Was not in burst.
                    if Criteria(i) % Criteria met, now in new burst.
                        INBURST = 1; % Update.
                        NUM_ = NUM_ + 1;
                        NUMBER = NUM_;
                        BL = 1;
                    else % Still not in burst, continue.
                        continue
                    end
                else % Was in burst.
                    if ~ Criteria(i) % Criteria no longer met.
                        INBURST = 0; % Update.
                        if BL<N % Erase if not big enough.
                            SpikeBurstNumber(SpikeBurstNumber==NUMBER) = -1;
                            NUM_ = NUM_ - 1;
                        end
                        NUMBER = -1;
                    elseif diff(obj.Spikes.Times([i-(N-1) i])) > ISI_N && BL >= N
                        % This conditional statement is necessary to split apart
                        % consecutive bursts that are not interspersed by a tonic spike
                        % (i.e. Criteria == 0). Occasionally in this case, the second
                        % burst has fewer than 'N' spikes and is therefore deleted in
                        % the above conditional statement (i.e. 'if BL<N').
                        %
                        % Skip this if at the start of a new burst (i.e. 'BL>=N'
                        % requirement).
                        %
                        NUM_ = NUM_ + 1; % New burst, update number.
                        NUMBER = NUM_;
                        BL = 1; % Reset burst length.
                    else % Criteria still met.
                        BL = BL + 1; % Update burst length.
                    end
                end
                SpikeBurstNumber(i) = NUMBER; % Assign a burst number to
                % each spike.
            end
            % %% Assign Burst information
            MaxBurstNumber = max(SpikeBurstNumber);
            Burst.T_start = zeros(1,MaxBurstNumber); % Burst start time [sec]
            Burst.T_end = zeros(1,MaxBurstNumber); % Burst end time [sec]
            for i = 1:MaxBurstNumber
                ID = find( SpikeBurstNumber==i );
                Burst.T_start(i) = obj.Spikes.Times(ID(1));
                Burst.T_end(i) = obj.Spikes.Times(ID(end));
            end
            run_time = toc;
            fprintf('Detected %i bursts in %.2f seconds using %0.2f minutes of spike data.\n', ...
                MaxBurstNumber,run_time,diff(obj.Spikes.Times([1 end]))/60);
        end
        
    end
    
    methods (Static)
        
        function ParentPath = getParentPath()
            ParentPath = fileparts(fileparts(which('MEArecording')));
        end
        
        function defaultParams = returnDefaultParams()
            % Fields to set true for the respective analyses that should be
            % performed
            defaultParams.Analyses.SingleCell = true;
            defaultParams.Analyses.Regularity = true;
            defaultParams.Analyses.Bursts = true;
            defaultParams.Analyses.Connectivity.CCG = false;
            defaultParams.Analyses.Connectivity.DDC = true;
            
            % Parameters for outlier detection 
            defaultParams.Outlier.Method = "median"; %Method as defined by Matlab's rmoutliers() function // set to [] to skip outlier detection
            defaultParams.Outlier.ThresholdFactor = 3; %Corresponding threshold
            
            % Parameters for the unit quality control (QC)
            defaultParams.QC.LSB = 6.2; %Factor by which to multiply voltage data
            defaultParams.QC.Amplitude = [20 300]; %[Minimum Maximum] amplitude
            defaultParams.QC.FiringRate = [0.1 10]; %[Minimum Maximum] firing rate
            defaultParams.QC.RPV = 0.02; %Refractory period violation ratio
            defaultParams.QC.Axon = 0.8; %Ratio (positive amplitude)/ (maximum amplitude) to be considered axonal
            defaultParams.QC.Noise = 14; % # of sign changes for the waveform derivative ("up-and-down" waveform)
            defaultParams.QC.NoiseCutout = [-1 1]; % [ms] before and after peak to be considered for noise detection
            
            % Parameters for the burst detection
            defaultParams.Bursts.N = 0.9; %0.0015; %Number of spikes in bursts (if < 1 used as a percentage of the mean network firing rate)
            defaultParams.Bursts.N_min = 65; %Minimum value for N if used as a ratio (N < 1)
            defaultParams.Bursts.ISI_N = 1.5; %Relates to minimum burst length [s]
            defaultParams.Bursts.Merge_t = 2; % Maximum ratio of IBI/BD to merge bursts to prevent oversplitting
            defaultParams.Bursts.Binning = 0.1;
            
            % Parameters for the burst detection
            defaultParams.Regularity.Binning = 0.1; %Binning to compute regularity features
            defaultParams.Regularity.N_peaks = 20; % Number of peaks for regularity fit
            
            % Parameters for CCG calculation
            defaultParams.CCG.BinSize = .001; %1ms
            defaultParams.CCG.Duration = 0.1; %100ms, corresponds to the length of one side of ccg
            defaultParams.CCG.Conv_w = .010/defaultParams.CCG.BinSize;  % 10ms window (gaussian convolution)
            defaultParams.CCG.Alpha = 0.001; %Threshold for significance testing
            
            % Parameters for DDC calculation
            defaultParams.DDC.BinSize = .001; %1ms
            defaultParams.DDC.Threshold = 0; %ReLU treshold
            
        end
    end
    
    methods 
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Analyses
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function unit_array = generateUnits(obj)
           [max_amplitudes, reference_electrode, norm_wf_matrix] = obj.generateWaveformMatrix();
           [firing_rates, unit_spike_times] = obj.calculateFiringRates();
           good_units = performUnitQC(obj, max_amplitudes, firing_rates, unit_spike_times, norm_wf_matrix);
           if sum(good_units) > 0
               good_amplitudes = max_amplitudes(good_units);
               good_unit_spike_times = unit_spike_times(good_units);
               good_wf_matrix = norm_wf_matrix(:,good_units);
               waveform_features = obj.inferWaveformFeatures(good_amplitudes, good_wf_matrix);
               unit_array = Unit();
               for u = 1:length(waveform_features)
                   unit_array(u) = Unit(obj, good_wf_matrix(:,u), reference_electrode(u), good_unit_spike_times{u}, waveform_features{u});
               end
               obj.Units = unit_array;
               obj.updateSpikeTimes();
               obj.aggregateSingleCellFeatures();
           else
               warning("No good units found")
           end
       end
       
       function getRegularity(obj)
           binned_network_activity = histcounts(obj.Spikes.Times,'BinLimits',[0 obj.RecordingInfo.Duration],'BinWidth',obj.Parameters.Regularity.Binning);
           norm_nw_activity = binned_network_activity/max(binned_network_activity);
           NFFT = length(norm_nw_activity);
           F = (0 : 1/NFFT : 1/2-1/NFFT)*(1/obj.Parameters.Regularity.Binning);
           TEMP = fft(norm_nw_activity,NFFT);
           TEMP(1) = 0;
           freq_domain = abs(TEMP(1:NFFT/2));
           [mag,freq_idx] = max(freq_domain);
           obj.NetworkFeatures.Regularity.NetworkRegularityFrequency = F(freq_idx);
           obj.NetworkFeatures.Regularity.NetworkRegularityMagnitude = mag;
           
           norm_freq_domain = freq_domain/max(freq_domain);
           l = [];
           p = [];
           for i = 1:obj.Parameters.Regularity.N_peaks
               [pks,locs,~,~] = findpeaks(norm_freq_domain,'NPeaks',1,'SortStr','descend');
               norm_freq_domain = norm_freq_domain(locs:end);
               l = [l; locs]; %#ok<AGROW>
               p = [p; pks]; %#ok<AGROW>
               if isempty(pks)
                   break
               end
           end
           log_p = log10(p)-min(log10(p)); % Make data positive to be able to fit
           try
               f = fit(cumsum(l),log_p,'exp1');
               obj.NetworkFeatures.Regularity.NetworkRegularityFit = f.b;
           catch
               obj.NetworkFeatures.Regularity.NetworkRegularityFit = NaN;
           end
       end
       
       function Burst = detectBursts(obj)
           if obj.Parameters.Bursts.N < 1
               N = ceil((length(obj.Spikes.Times) / obj.RecordingInfo.Duration) * obj.Parameters.Bursts.N);
           else
               N = obj.Parameters.Bursts.N;
           end
           
           N = max([obj.Parameters.Bursts.N_min N]); %Minimum of 65 spikes per burst // empirical minimum
           %             ISI_N = (N > 65) * 0.5 + 1; %Reduce ISI_N if N = 65
           ISI_N = obj.Parameters.Bursts.ISI_N;
           IBI_merge = obj.Parameters.Bursts.Merge_t * ISI_N;
           Burst = detectBurstsISIN(obj, N, ISI_N);
           N_bursts = length(Burst.T_start);
           
           %Merge bursts with too short IBIs (< IBI_merge)
           short_burst_idx = find(Burst.T_start(2:end) - Burst.T_end(1:end-1)> IBI_merge);
           Burst.T_start = Burst.T_start([1 short_burst_idx+1]);
           Burst.T_end = Burst.T_end([short_burst_idx length(Burst.T_end)]);
           fprintf('Merged %i bursts\n', N_bursts - length(Burst.T_start))
           obj.Bursts = Burst;
       end
       
       function getBurstStatistics(obj)
          if isempty(obj.Bursts)
              detectBursts(obj);
          end
          
          IBIs = obj.Bursts.T_start(2:end) - obj.Bursts.T_end(1:end-1);
          BDs = obj.Bursts.T_end - obj.Bursts.T_start;
          
          [RiseTime, FallTime, RiseRate, FallRate] = deal(nan(1,length(BDs)));
          Fs = round(1/obj.Parameters.Bursts.Binning);
          inburst_spikes = 0;
          for iburst = 1:length(BDs)
              burst_activity = obj.Spikes.Times(obj.Spikes.Times <= obj.Bursts.T_end(iburst) & obj.Spikes.Times >= obj.Bursts.T_start(iburst));
              binned_activity = histcounts(burst_activity,'BinWidth',obj.Parameters.Bursts.Binning);
              [max_binned, imax_binned] = max(binned_activity);
              norm_activity = binned_activity/max_binned;
              pre_lower_level = min(norm_activity(1:imax_binned));
              post_lower_level = min(norm_activity(imax_binned:end));
              RiseTime(iburst) = risetime([ones(1,Fs)*pre_lower_level norm_activity(1:imax_binned) ones(1,Fs)],Fs);
              FallTime(iburst) = falltime([ones(1,Fs) norm_activity(imax_binned:end) ones(1,Fs)*post_lower_level],Fs);
              RiseRate(iburst) = 0.8/RiseTime(iburst);
              FallRate(iburst) = 0.8/FallTime(iburst);
              inburst_spikes = inburst_spikes + length(burst_activity);
          end
          
          cellfun_input = {IBIs, BDs, RiseTime, FallTime, RiseRate, FallRate};
          
          % Perform outlier removal if method is specified
          if ~isempty(obj.Parameters.Outlier.Method)
              cellfun_input = cellfun(@(x) rmoutliers(x,obj.Parameters.Outlier.Method,'ThresholdFactor',obj.Parameters.Outlier.ThresholdFactor),...
                  cellfun_input,'un',0);
          end
          feature_means = cellfun(@mean, cellfun_input,'un',0);
          [obj.NetworkFeatures.Burst.MeanInterBurstInterval,...
              obj.NetworkFeatures.Burst.MeanBurstDuration,...
              obj.NetworkFeatures.Burst.MeanRiseTime,...
              obj.NetworkFeatures.Burst.MeanFallTime,...
              obj.NetworkFeatures.Burst.MeanRiseVelocity,...
              obj.NetworkFeatures.Burst.MeanDecayVelocity] = feature_means{:};
          obj.NetworkFeatures.Burst.VarianceBurstDuration = var(BDs);
          obj.NetworkFeatures.Burst.VarianceInterBurstInterval = var(IBIs);
          obj.NetworkFeatures.Burst.IntraBurstFiringRate = inburst_spikes/obj.RecordingInfo.Duration;
          obj.NetworkFeatures.Burst.InterBurstFiringRate = (length(obj.Spikes.Times) - inburst_spikes)/obj.RecordingInfo.Duration;
       end
       
       function [ccg,t] = calculateCCGs(obj)
           if isempty(which('CCGHeart'))
               cd(fullfile(obj.getParentPath(),'Functions'))
               mex('CCGHeart.c')
           end
           [ccg,t] = CCG(obj.Spikes.Times,obj.Spikes.Units,'binSize',obj.Parameters.CCG.BinSize,'duration',obj.Parameters.CCG.Duration,'Fs',1/obj.RecordingInfo.SamplingRate);
           disp('Finished CCG calculations')
       end
       
       function inferConnectivityCCG(obj)
           [ccgR1,tR] = obj.calculateCCGs();
           binSize = obj.Parameters.CCG.BinSize; %.5ms
           conv_w = obj.Parameters.CCG.Conv_w;  % 10ms window
           alpha = obj.Parameters.CCG.Alpha; %high frequency cut off, must be .001 for causal p-value matrix
%            Fs = 1/obj.RecordingInfo.SamplingRate;
           
           nCel=size(ccgR1,2);
           
           % Create CCGs (including autoCG) for all cells
%            [ccgR1,tR] = CCG(sorted_spk_vec,sorted_spk_idx,'binSize',binSize,'duration',duration,'Fs',Fs);
           
           ccgR = nan(size(ccgR1,1),nCel,nCel);
           ccgR(:,1:size(ccgR1,2),1:size(ccgR1,2)) = ccgR1;
           
           % get  CI for each CCG
           Pval=nan(length(tR),nCel,nCel);
           Pred=zeros(length(tR),nCel,nCel);
           Bounds=zeros(size(ccgR,1),nCel,nCel);
           sig_con = [];
           sig_con_inh = [];
           
           for refcellID=1:max(nCel)
               for cell2ID=1:max(nCel)
                   
                   if(refcellID==cell2ID)
                       continue;
                   end
                   
                   cch=ccgR(:,refcellID,cell2ID);			% extract corresponding cross-correlation histogram vector
                   
                   [pvals,pred,~]=bz_cch_conv(cch,conv_w);
                   % Store predicted values and pvalues for subsequent plotting
                   Pred(:,refcellID,cell2ID)=pred;
                   Pval(:,refcellID,cell2ID)=pvals(:);
                   Pred(:,cell2ID,refcellID)=flipud(pred(:));
                   Pval(:,cell2ID,refcellID)=flipud(pvals(:));
                   
                   % Calculate upper and lower limits with bonferonni correction
                   % monosynaptic connection will be +/- 4 ms
                   
                   nBonf = round(.005/binSize)*2;
                   hiBound=poissinv(1-alpha/nBonf,pred);
                   loBound=poissinv(alpha/nBonf, pred);
                   Bounds(:,refcellID,cell2ID,1)=hiBound;
                   Bounds(:,refcellID,cell2ID,2)=loBound;
                   
                   Bounds(:,cell2ID,refcellID,1)=flipud(hiBound(:));
                   Bounds(:,cell2ID,refcellID,2)=flipud(loBound(:));
                   
                   sig = cch>hiBound;
                   sig_inh= cch < loBound;
                   
                   % Find if significant periods falls in monosynaptic window +/- 4ms
                   prebins = round(length(cch)/2 - .0032/binSize):round(length(cch)/2);
                   postbins = round(length(cch)/2 + .0008/binSize):round(length(cch)/2 + .004/binSize);
                   cchud  = flipud(cch);
                   sigud  = flipud(sig);
                   sigud_inh=flipud(sig_inh);
                   
                   sigpost=max(cch(postbins))>poissinv(1-alpha,max(cch(prebins)));
                   sigpre=max(cchud(postbins))>poissinv(1-alpha,max(cchud(prebins)));
                   
                   sigpost_inh=min(cch(postbins))<poissinv(alpha,mean(cch(prebins)));
                   sigpre_inh=min(cchud(postbins))<poissinv(alpha,mean(cchud(prebins)));
                   %check which is bigger
                   if (any(sigud(prebins)) && sigpre)
                       %test if causal is bigger than anti causal
                       sig_con = [sig_con;cell2ID refcellID];
                   end
                   
                   if (any(sig(postbins)) && sigpost)
                       sig_con = [sig_con;refcellID cell2ID];
                   end
                   
                   if (any(sigud_inh(prebins)) && sigpre_inh)
                       %test if causal is bigger than anti causal
                       sig_con_inh = [sig_con_inh;cell2ID refcellID];
                   end
                   if (any(sig_inh(postbins)) && sigpost_inh)
                       sig_con_inh = [sig_con_inh;refcellID cell2ID];
                   end
               end
               
           end
           
           sig_con=unique(sig_con,'rows');
           sig_con_inh=unique(sig_con_inh, 'rows');
           
           if(~isempty(sig_con))
               ccg_vec=[];
               for jj=1:size(sig_con,1)
                   ccg_vec=[ccg_vec ccgR(:,sig_con(jj,1),sig_con(jj,2))];
               end
           else
               ccg_vec=[];
           end
           
           if(~isempty(sig_con_inh))
               ccg_vec_inh=[];
               for jj=1:size(sig_con_inh,1)
                   ccg_vec_inh=[ccg_vec_inh ccgR(:,sig_con_inh(jj,1),sig_con_inh(jj,2))];
               end
           else
               ccg_vec_inh=[];
           end
           fprintf('Found %i excitatory and %i inhibitory connections\n', size(sig_con,1), size(sig_con_inh,1))
           con_mat = zeros(size(ccgR,[2,3]));
           for e = 1:size(sig_con,1)
               con_mat(sig_con(e,1),sig_con(e,2)) = 1;
           end
           for i = 1:size(sig_con_inh,1)
               con_mat(sig_con_inh(i,1),sig_con_inh(i,2)) = -1;
           end
           obj.Connectivity.CCGResults.ExcitatoryConnection = sig_con;
           obj.Connectivity.CCGResults.InhibitoryConnection = sig_con_inh;
           obj.Connectivity.CCGResults.ConnectivityMatrix = con_mat;
           obj.Connectivity.CCGResults.CCGs = ccgR;
       end
       
       function inferConnectivityDDC(obj)%[Cov,precision,B,dCov]
           %{
            INPUT:
            V_obs: time points x variables
            thres: Relu offset
            TR: sampling interval (seconds)
            OUTPUT:
            B: <ReLu(x),x>
            dCov: <dx/dt,x>
            Version2: change the implementation of ReLu, get rid of the translation
           %}
           TR = 1;
           V_obs = zeros(obj.RecordingInfo.Duration * round(1/obj.Parameters.DDC.BinSize),length(obj.Units));
           spiketimes_ms = round(obj.Spikes.Times * 1000);
           V_obs(sub2ind(size(V_obs), spiketimes_ms, obj.Spikes.Units)) = 1;
           V_obs = smoothdata(V_obs,'gaussian',3);
           [~,N] = size(V_obs);
           % 	Cov = cov(V_obs);
           % 	precision = inv(Cov);
           Fx = V_obs-obj.Parameters.DDC.Threshold; Fx(Fx<0) = 0;
           % Fx = V_obs; Fx(Fx<thres) = 0;
           tmp = cov([Fx V_obs]);
           B = tmp(1:N,N+1:end);
           dV = (-1/2*V_obs(1:end-2,:) + 1/2*V_obs(3:end,:))/TR; % (z(t+1)-z(t-1))/2
           dV = [mean(dV);dV;mean(dV)]; % substitute the first and last row with mean
           tmp = cov([dV V_obs]);
           dCov = tmp(1:N,N+1:end);
           obj.Connectivity.DDC = dCov * inv(B);
       end
       
       function feat_table = concatenateFeatures(~,feat_struct, feat_names) %Usable for unit and network features
           feat_array = squeeze(struct2cell(feat_struct));
           if ~isequal(feat_names,"all")
               feat_array = feat_array(contains(fieldnames(feat_struct), feat_names,'IgnoreCase',true),:);
           end
           feat_tables = cellfun(@struct2table, feat_array,'un',0);
           if isstruct(feat_tables{1}(1,1).Variables) %For unit features
               feat_tables = arrayfun(@(x) struct2table(feat_tables{1}(1,x).Variables),1:size(feat_tables{1},2),'un',0)';
           end
           feat_table = arrayfun(@(x) [feat_tables{:,x}],1:size(feat_tables,2),'un',0);
           feat_table = vertcat(feat_table{:});
       end
       
       function feat_table = getFeatureTable(obj, unit_features, network_features)
           arguments
               obj MEArecording
               unit_features string = "all"
               network_features string = "all"
           end
           unit_table = concatenateFeatures(obj,[obj.UnitFeatures], unit_features);
           network_table = concatenateFeatures(obj,[obj.NetworkFeatures], network_features);
           feat_table = [unit_table network_table];
       end
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Plots
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function PlotNetworkScatter(obj)
          figure('Color','w');
          plot(obj.Spikes.Times,obj.Spikes.Units,'k.')
          xlabel('Time [s]')
          ylabel('Unit ID')
          axis tight
          box off
       end
       
       function PlotNetworkScatterHistogram(obj)
           N_bins = [round(obj.RecordingInfo.Duration / obj.Parameters.Bursts.Binning); length(obj.Units)];
           figure('Color','w');
           scatterhistogram(obj.Spikes.Times,obj.Spikes.Units,'MarkerSize',1,'Color','k','HistogramDisplayStyle','stairs','NumBins', N_bins);
           
       end
       
       function PlotBurstCheck(obj)
           figure('Color','w');
           histogram(obj.Spikes.Times,'BinWidth',obj.Parameters.Bursts.Binning,'DisplayStyle','stairs')
           if ~isempty(obj.Bursts)
               hold on
               xline(obj.Bursts.T_start,'g--','LineWidth',1,'Label','Start','LabelHorizontalAlignment','center')
               xline(obj.Bursts.T_end,'r--','LineWidth',1,'Label','End','LabelHorizontalAlignment','center')
           else
               warning('No burst data found. Run obj.detectBursts() first')
           end
           box off
       end
       
       function PlotCCG(obj,ccg_idx1,ccg_idx2)
           f = figure('Color','w');
           a1 = axes(f);
           if ccg_idx1 ~= ccg_idx2
               [ccg,t] = CCG({obj.Units(ccg_idx1).SpikeTimes,obj.Units(ccg_idx2).SpikeTimes},[],'binSize',obj.Parameters.CCG.binSize,...
                   'duration',obj.Parameters.CCG.duration,'Fs',1/obj.RecordingInfo.SamplingRate);
               bar(a1,t,ccg(:,1,2),'k','BarWidth',1)
           else
               [ccg,t] = CCG(spk_data.spike_times(ccg_idx1),[],'binSize',obj.Parameters.CCG.binSize,...
                   'duration',obj.Parameters.CCG.duration,'Fs',1/obj.RecordingInfo.SamplingRate);
               bar(a1,t,ccg,'k','BarWidth',1)
           end
           set(gca,'FontSize',7)
                      
           xlabel('Interval (ms)')
           xline(0,'w');
           % xlabel('Interval [ms]')
           a1.YAxis.Visible = 'off';
           % xticks([-25])
           box off
           axis tight
           xlims = xlim;
           xticks([xlims(1) 0 xlims(2)])
           x_max = obj.Parameters.CCG.duration/2*1000;
           xticklabels([-x_max 0 x_max])
       end
       
       function indsort = communityPlot(con_mat,flag)
           %flag indicates if squares are drawnaround communities
           % M     = community_louvain(con_mat,0.0001,[],'negative_asym');
           M = modularity_dir(con_mat,100);
           [X,Y,indsort] = grid_communities(M);
           imagesc(con_mat(indsort,indsort))
           if flag
               hold on
               plot(X,Y,'k','linewidth',1);
           end
           cm = othercolor('RdBu5',3);
           colormap(flipud(cm))
       end
    end
end