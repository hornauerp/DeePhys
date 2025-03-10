classdef MEArecording < handle
    
    %%%%%%%TODO%%%%%%%
    % 
    %%%%%%%%%%%%%%%%%%
    
    properties
        Metadata
        RecordingInfo
        Parameters
        Spikes
        Bursts
        Units
        NetworkFeatures
        UnitFeatures
        ClusteredFeatures
        Connectivity
        NumUnitClusters % Number of unit clusters
    end

    properties (Dependent)
        ParentRecordingPath
        FullACGs
        ACGs
    end

    methods (Hidden)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Object generation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        function mearec = MEArecording(metadata, parameters)
            arguments
                metadata (1,1) struct = struct() %Metadata to be parsed // metadata.InputPath is required
                parameters (1,1) struct = struct()
            end
            
            if nargin > 0
                addpath(genpath(mearec.getParentPath()));
                mearec.parseMetadata(metadata);
                mearec.parseRecordingInfo();
                mearec.parseParameters(parameters);
                mearec.performAnalyses();
                mearec.saveObject();
            end
        end

                
        function parseMetadata(obj, metadata)
            arguments
               obj MEArecording
               metadata struct
            end
            
            if isempty(fieldnames(metadata)) %Check if any metadata is provided
                error("No metadata provided, cannot continue");
            elseif ~isfield(metadata,"InputPath")
                error("Metadata does not provide InputPath, cannot continue");
            elseif isfield(metadata,"LookupPath") && ~isempty(metadata.LookupPath)
                if ~isfile(metadata.LookupPath)
                    error("LookUpPath does not link to a valid file path")
                end
                obj.Metadata = obj.retrieveMetadata(metadata);
                disp('Import metadata from file')
            else
                obj.Metadata = metadata;
                warning("Used metadata provided, no import")
            end
            
            if isfield(obj.Metadata,"RecordingDate") && ~isempty(obj.Metadata.RecordingDate) && isfield(obj.Metadata,"PlatingDate") && ~isempty(obj.Metadata.PlatingDate)
               obj.Metadata.DIV = days(datetime(obj.Metadata.RecordingDate,"InputFormat","yyMMdd") - datetime(num2str(obj.Metadata.PlatingDate),"InputFormat","yyMMdd"));
            end
        end
        
        function parseRecordingInfo(obj)
            obj.RecordingInfo.ElectrodeCoordinates = obj.getElectrodeCoordinates();
            obj.RecordingInfo.SamplingRate = obj.getSamplingRate();
            [obj.Spikes.Times, obj.Spikes.Units, obj.Spikes.hasSpikes] = obj.getSpikeTimes();
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
        
        function [spike_times, spike_units, hasSpikes] = getSpikeTimes(obj)
            spike_times = double(readNPY(fullfile(obj.Metadata.InputPath, "spike_times.npy"))) / obj.RecordingInfo.SamplingRate;
            spike_units = double(readNPY(fullfile(obj.Metadata.InputPath, "spike_templates.npy")))+1;
            
            % Find unit ids that do not have any spikes
            max_ids = 1:max(spike_units);
            hasSpikes = ismember(max_ids,spike_units);
        end
        
        function template_matrix = getTemplateMatrix(obj)
            template_matrix = readNPY(fullfile(obj.Metadata.InputPath, "templates.npy"));
%             template_matrix = template_matrix(obj.Spikes.hasSpikes,:,:);
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
                for pf = 1:length(parameter_fields)
                    if ismember(parameter_fields(pf),fieldnames(obj.Parameters))
                        subfields = fieldnames(parameters.(parameter_fields{pf}));
                        for sf = 1:length(subfields)
                            if ismember(subfields{sf},fieldnames(obj.Parameters.(parameter_fields{pf})))
                                obj.Parameters.(parameter_fields{pf}).(subfields{sf}) = parameters.(parameter_fields{pf}).(subfields{sf});
                            else
                                warning("Unrecognized parameter field: %s.%s",parameter_fields{pf},subfields{sf})
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
                if isempty(obj.Units)
                    return
                end
                obj.UnitFeatures = obj.aggregateSingleCellFeatures(obj.Units);
            end
            
            if obj.Parameters.Analyses.Regularity
                obj.getRegularity();
            end
            
            if obj.Parameters.Analyses.Catch22
                obj.NetworkFeatures.Catch22 = obj.run_catch_22();
            end
            
            if obj.Parameters.Analyses.Bursts
                obj.getBurstStatistics();
            end
            
            if ~isempty(obj.Parameters.Analyses.Connectivity)
                obj.inferConnectivity(obj.Parameters.Analyses.Connectivity);
                if ~isempty(fields(obj.Connectivity))
                    obj.inferGraphFeatures();
                end
            end
        end
        
        function saveObject(obj)
            if obj.Parameters.Save.Flag
                if isempty(obj.Parameters.Save.Path)
                    save_path = fullfile(obj.Metadata.InputPath, 'MEArecording.mat');
                else
                    save_path = fullfile(obj.Parameters.Save.Path, 'MEArecording.mat');
                end
                
                if exist(save_path,'file') && ~obj.Parameters.Save.Overwrite
                    return
                else
                    save(save_path, "obj")
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Quality control
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [max_amplitudes, reference_electrode, norm_wf_matrix] = generateWaveformMatrix(obj)
            template_matrix = obj.getTemplateMatrix();
            template_matrix = template_matrix(:,sum(template_matrix,[1,3],'omitnan')~=0,:); %Remove all zero paddings
            [max_amplitudes, max_idx] = min(template_matrix,[],[2,3],'linear');
            max_amplitudes = abs(max_amplitudes * obj.Parameters.QC.LSB);
            [~,~,reference_electrode] = ind2sub(size(template_matrix),max_idx);
            peak_wf = arrayfun(@(x) template_matrix(x,:,reference_electrode(x)),1:length(reference_electrode),'un',0);
            peak_wf_matrix = cat(1,peak_wf{:})';
            norm_wf_matrix = peak_wf_matrix./max(abs(peak_wf_matrix));
        end
        
        function [firing_rates, unit_spike_times] = calculateFiringRates(obj, N_units)
            [spike_times, spike_units] = obj.getSpikeTimes();
            unit_spike_times = arrayfun(@(x) spike_times(spike_units == x), 1:N_units,'un',0);
            firing_rates = cellfun(@length, unit_spike_times)/obj.RecordingInfo.Duration;
        end
        
        function good_units = performUnitQC(obj, max_amplitudes, firing_rates, unit_spike_times, norm_wf_matrix)
            bad_power = checkPower(obj, norm_wf_matrix);
            bad_amplitude = checkAmplitude(obj, max_amplitudes);
            bad_firing_rate = checkFiringRate(obj,firing_rates);
            bad_rpv = checkRefractoryPeriodViolations(obj, unit_spike_times);
            [is_axonal, is_noisy] = checkWaveform(obj, norm_wf_matrix);
            good_units = ~(bad_power | bad_amplitude | bad_firing_rate | bad_rpv | is_axonal | is_noisy);
            fprintf('Found %i/%i good units\n',sum(good_units),length(good_units))
        end
        
        function bad_power = checkPower(obj, norm_wf_matrix)
            
            n = size(norm_wf_matrix,1);
            y = fft(norm_wf_matrix);
            power = abs(y).^2/n;
            max_p = max(power)';
            bad_power = max_p > obj.Parameters.QC.PowerCutoff;
            fprintf('Found %i units with bad power\n',sum(bad_power))
        end
        
        function bad_amplitude = checkAmplitude(obj, max_amplitudes)
            if ~isempty(obj.Parameters.QC.Amplitude) && all(~isnan(obj.Parameters.QC.Amplitude))
                bad_amplitude = max_amplitudes < obj.Parameters.QC.Amplitude(1) | max_amplitudes > obj.Parameters.QC.Amplitude(2) | isnan(max_amplitudes);
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
            if ~isempty(obj.Parameters.QC.Axon) && ~isnan(obj.Parameters.QC.Axon)
                amplitude_ratio = -(max(norm_wf_matrix)./min(norm_wf_matrix))';
                %Fit two gaussians and try to find cutoff
                is_axonal = amplitude_ratio > obj.Parameters.QC.Axon;
            else
                is_axonal = zeros(size(norm_wf_matrix,2),1);
            end
            fprintf('Identified %i units as axonal\n',sum(is_axonal))
            
            if ~isempty(obj.Parameters.QC.Noise) && ~isnan(obj.Parameters.QC.Noise)
                [~,peak_amp_idx] = min(norm_wf_matrix);
                avg_peak_idx = median(peak_amp_idx(peak_amp_idx >1));
                cutout_idx = round(avg_peak_idx) + obj.Parameters.QC.NoiseCutout(1) * 10:...
                    round(avg_peak_idx) + obj.Parameters.QC.NoiseCutout(2) * 10;
                bad_peak = peak_amp_idx < cutout_idx(1) | peak_amp_idx > cutout_idx(end);
                noise_indicator = sum(abs(diff(diff(norm_wf_matrix(cutout_idx,:))>0)))'; %Check for non-monotomies in waveform shapes as noise indicators
                
                is_noisy = noise_indicator > obj.Parameters.QC.Noise;
                is_noisy = is_noisy | bad_peak';
            else
                is_noisy = zeros(size(norm_wf_matrix,2),1);
            end
            fprintf('Identified %i units as noise\n',sum(is_noisy))
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analyses
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function waveform_features = inferWaveformFeatures(obj, max_amplitudes, norm_wf_matrix)
            interpolation_factor = 10;
            ms_conversion = 10 * interpolation_factor; %Need to find a way to automate that, maybe 10/ms is standard for Ph?(obj.RecordingInfo.SamplingRate / 1000) * interpolation_factor;
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
                zc_pre_trough = max([zc_pre_trough 2]);
                unit_features.AUC_peak_1 = trapz(interp_wf_matrix(zc(zc_pre_trough - 1):zc(zc_pre_trough)));
                unit_features.AUC_trough = abs(trapz(interp_wf_matrix(zc(zc_pre_trough):zc(zc_pre_trough + 1))));
                if length(zc(zc_pre_trough:end)) > 3 %Ensure that there is another zero crossing at the end
                    unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):zc(zc_pre_trough + 2)));
                else %Otherwise we integrate until the end of the waveform
                    unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):size(interp_wf_matrix,1)));
                end
                %We pad the signal with min/max values to ensure reliable slewrate
                %calculations
                rise_cutout = interp_wf_matrix(unit_trough_idx(u):peak_2_idx(u),u);
                padded_rise = [ones(100,1)*rise_cutout(1); rise_cutout; ones(100,1)*rise_cutout(end)];
                unit_features.Rise = mean(slewrate(padded_rise,interpolation_factor));
                decay_cutout = interp_wf_matrix(peak_2_idx(u):end,u);
                decay_cutout = decay_cutout(1:find(decay_cutout==min(decay_cutout)));
                padded_decay = [ones(100,1)*decay_cutout(1); decay_cutout; ones(100,1)*decay_cutout(end)];
                unit_features.Decay = mean(slewrate(padded_decay,interpolation_factor));
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
        
        function aggregated_struct = aggregateSingleCellFeatures(obj, unit_array, features)
            arguments
                obj MEArecording
                unit_array Unit = [obj.Units]
                features string = ["ActivityFeatures","WaveformFeatures","RegularityFeatures","Catch22"];%,"GraphFeatures"]
            end
            
            for f = 1:length(features)
                if length(unit_array) > 1
                    feature_table = vertcat(unit_array.(features(f)));
                    %
                    if ~isempty(obj.Parameters.Outlier.Method)
                        feature_means = mean(rmoutliers(feature_table.Variables, obj.Parameters.Outlier.Method,'ThresholdFactor',obj.Parameters.Outlier.ThresholdFactor),1,"omitnan");
                    else
                        feature_means = mean(feature_table.Variables,1,"omitnan");
                    end
                    
                    aggregated_struct.(features(f)) = array2table(feature_means,"VariableNames",feature_table.Properties.VariableNames);
                else
                    aggregated_struct.(features(f)) = unit_array.(features(f));
                end
            end
        end
        
        function calculateClusterSingleCellFeatures(obj)
            for c = 1:obj.NumUnitClusters
                unit_array = obj.Units([obj.Units.ClusterID] == c);
                if isempty(unit_array) %Add empty clusters to have consistent feature matrix sizes
                    fnames = fieldnames(obj.UnitFeatures);
                    for f = 1:length(fnames)
                    	empty_struct.(fnames{f}) = obj.UnitFeatures.(fnames{f});
                        empty_struct.(fnames{f}).Variables = zeros(size(empty_struct.(fnames{f}).Variables));
                    end
                    obj.ClusteredFeatures{c} = empty_struct;
                else
                    obj.ClusteredFeatures{c} = obj.aggregateSingleCellFeatures(unit_array);
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
            if N == 0
                Burst.T_start = [];
                Burst.T_end = [];
            else
                tic;
                [spike_times, ~] = obj.remove_tonic_units();
                dT = zeros(N,length(spike_times))+inf;
                for j = 0:N-1
                    dT(j+1,N:length(spike_times)-(N-1)) = spike_times( (N:end-(N-1))+j ) - ...
                        spike_times( (1:end-(N-1)*2)+j );
                end
                Criteria = zeros(size(spike_times)); % Initialize to zero
                Criteria( min(dT)<=ISI_N ) = 1; % Spike passes condition if it is
                % included in a set of N spikes
                % with ISI_N <= threshold.
                % %% Assign burst numbers to each spike
                SpikeBurstNumber = zeros(size(spike_times)) - 1; % Initialize to '-1'
                INBURST = 0; % In a burst (1) or not (0)
                NUM_ = 0; % Burst Number iterator
                NUMBER = -1; % Burst Number assigned
                BL = 0; % Burst Length
                for i = N:length(spike_times)
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
                        elseif diff(spike_times([i-(N-1) i])) > ISI_N && BL >= N
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
                    Burst.T_start(i) = spike_times(ID(1));
                    Burst.T_end(i) = spike_times(ID(end));
                end
                run_time = toc;
                fprintf('Detected %i bursts in %.2f seconds using %0.2f minutes of spike data.\n', ...
                    MaxBurstNumber,run_time,diff(spike_times([1 end]))/60);
            end
        end
        
        function empirical_sttc = empiricalSTTCmatrix(obj)
            empirical_sttc = zeros(length(obj.Units));
            for i = 1:length(obj.Units)-1
                N1v = length(obj.Units(i).SpikeTimes);
                for j = (i+1):length(obj.Units)
                        N2v = length(obj.Units(j).SpikeTimes);
                        empirical_sttc(i,j) = sttc(N1v, N2v, obj.Parameters.STTC.MaxLag, [0 obj.RecordingInfo.Duration], obj.Units(i).SpikeTimes, obj.Units(j).SpikeTimes);
                end
            end
            empirical_sttc = empirical_sttc+triu(empirical_sttc,-1).'; %Make matrix symmetric
        end
        
        function surrogate_sttc = surrogateSTTCmatrix(obj)
            silent_units = find(arrayfun(@(x) isempty(x.SpikeTimes), obj.Units));
            spks_data = convert2asdf2(obj.Spikes.Times, obj.Spikes.Units, obj.RecordingInfo.SamplingRate);
            
            for su = silent_units %Add empty units in case of concatenated recordings
                spks_data.raster = [spks_data.raster(1:su-1); cell(1,1); spks_data.raster(su:end)];
            end
            spks_data.nchannels = length(spks_data.raster);
            
            spks_rand = randomizeasdf2(spks_data, obj.Parameters.STTC.SurrogateMethod);
            surrogate_sttc = zeros(spks_data.nchannels);
            for i = 1:(spks_data.nchannels - 1)
                spike_times_1 = double(spks_rand.raster{i})/obj.RecordingInfo.SamplingRate;
                N1v = length(spike_times_1);
                for j = (i + 1):spks_data.nchannels
                    spike_times_2 = double(spks_rand.raster{j})/obj.RecordingInfo.SamplingRate;
                    N2v = length(spike_times_2);
                    surrogate_sttc(i,j) = sttc(N1v, N2v, obj.Parameters.STTC.MaxLag, [0 obj.RecordingInfo.Duration], spike_times_1, spike_times_2);
                end
            end
            surrogate_sttc = surrogate_sttc+triu(surrogate_sttc,-1).'; %Make matrix symmetric
        end
    end
    
    methods (Static)
        
        function ParentPath = getParentPath() %% Does not work when parallelized (e.g. parfor)
            ParentPath = fileparts(fileparts(which('MEArecording')));
        end
        
        function metadata_struct = retrieveMetadata(metadata)
            arguments
                metadata struct
                % Needs to contain indices for the path part that represents (in that order) 
                % metadata.PathPartIndex = [recording_date, plate_id, well_id]
            end
            
            metadata_struct = metadata;
            path_parts = strsplit(metadata.InputPath,"/");
            if ~isfield(metadata_struct,'RecordingDate')
                recording_date = path_parts(metadata.PathPartIndex(1));
                assert(~isnan(str2double(recording_date)),'Wrong PathPartIndex(1) / RecordingDate') %Check if recording date is a number
                metadata_struct.RecordingDate = recording_date;
            end
            
            if ~isfield(metadata_struct,'PlateID')
                plate_id = path_parts(metadata.PathPartIndex(2));
                plate_number = char(plate_id); 
                plate_number = plate_number(end-2:end);
                assert(~isnan(str2double(plate_number)),'Wrong PathPartIndex(2) / PlateID') %Check if recording date is a number
            else
                plate_id = metadata.PlateID;
            end
            
            well_id = char(path_parts(metadata.PathPartIndex(3)));
            well_id = well_id(end-1:end);
            assert(~isnan(str2double(well_id)),'Wrong PathPartIndex(3) / WellID')
            chip_id = plate_id + "_" + well_id;
            metadata_struct.PlateID = plate_id;
            metadata_struct.WellID = well_id;
            metadata_struct.ChipID = chip_id;

            info_table = readtable(metadata.LookupPath,"ReadVariableNames",true);
            for ids = 1:size(info_table,1)
                chip_ids = strtrim(string(strsplit(info_table.Chip_IDs{ids},",")));
                if contains(chip_id,chip_ids)
                    metadata_vars = string(info_table.Properties.VariableNames); 
                    metadata_vars = metadata_vars(~contains(metadata_vars,"Chip_IDs"));%Exclude Chip_IDs
                    for m = metadata_vars
                        try
                            metadata_struct.(m) = string(info_table.(m){ids});
                        catch
                            metadata_struct.(m) = info_table.(m)(ids);
                        end
                    end
                    break
                end
            end
            
        end
        
        function defaultParams = returnDefaultParams()
            % Fields to set true for the respective analyses that should be
            % performed
            defaultParams.Analyses.SingleCell = true;
            defaultParams.Analyses.Regularity = true;
            defaultParams.Analyses.Bursts = true;
            defaultParams.Analyses.Connectivity = ["CCG","STTC"]; %"CCG","STTC", "DDC" are implemented // can take several methods as input
            defaultParams.Analyses.Catch22 = true;
            
            % Parameters for outlier detection 
            defaultParams.Outlier.Method = "median"; %Method as defined by Matlab's rmoutliers() function // set to [] to skip outlier detection
            defaultParams.Outlier.ThresholdFactor = 5; %Corresponding threshold
            
            % Parameters for the unit quality control (QC)
            defaultParams.QC.LSB = 6.2; %Factor by which to multiply voltage data
            defaultParams.QC.Amplitude = [20 1000]; %[Minimum Maximum] amplitude
            defaultParams.QC.FiringRate = [0.01 10]; %[Minimum Maximum] firing rate
            defaultParams.QC.RPV = 0.02; %Refractory period violation ratio
            defaultParams.QC.Axon = 0.8; %Ratio (positive amplitude)/ (maximum amplitude) to be considered axonal
            defaultParams.QC.Noise = 14; % # of sign changes for the waveform derivative ("up-and-down" waveform)
            defaultParams.QC.NoiseCutout = [-1 1]; % [ms] before and after peak to be considered for noise detection
            defaultParams.QC.PowerCutoff = 1.3; % Maximum allowed power of the fourier transform of the reference waveform (excludes oscillating noise units)
            defaultParams.QC.N_Units = 10; % Minimum number of units that need to pass the QC to continue with the analysis
            defaultParams.QC.GoodUnits = []; % Units to keep if some previous analyses already determined good units IDs (e.g. manual curation, KS good units)
                                             % Skips the actual QC if not empty
            
            % Parameters for the burst detection
            defaultParams.Bursts.MergeFactor = 0.5; %Bursts with an IBI of MergeFactor * best_ISI_N will be merged
            defaultParams.Bursts.Binning = 0.1; %in [s] 0.1 = 100ms
            
            % Parameters for the burst detection
            defaultParams.Regularity.Binning = 0.1; %Binning to compute regularity features
            defaultParams.Regularity.N_peaks = 20; % Number of peaks for regularity fit
            
            % Parameters for CCG calculation
            defaultParams.CCG.BinSize = .001; %1ms
            defaultParams.CCG.Duration = 0.1; %100ms, corresponds to the length of one side of ccg
            defaultParams.CCG.Conv_w = .010/defaultParams.CCG.BinSize;  % 10ms window (gaussian convolution)
            defaultParams.CCG.Alpha = 0.001; %Threshold for significance testing
            
            % Parameters for GLM calculation
            defaultParams.GLM.binsize = 1;
            defaultParams.GLM.interval = 100;
            defaultParams.GLM.tau0 = 0.8;
            defaultParams.GLM.eta_w = 5;
            defaultParams.GLM.eta_tau = 20;
            defaultParams.GLM.eta_dt_coeff = 2;
            defaultParams.GLM.min_spikes = 100; 
            defaultParams.GLM.threshold = 5.09;
            
            % Parameters for DDC calculation
            defaultParams.DDC.BinSize = .001; %in seconds //1ms
            defaultParams.DDC.Threshold = 0; %ReLU treshold
            
            % Parameters for STTC calculation
            defaultParams.STTC.MaxLag = .05; %in seconds
            defaultParams.STTC.SurrogateMethod = 'jitter';
            defaultParams.STTC.N_Surrogates = 100;
            defaultParams.STTC.Percentile = 95;
            
            % Parameters for catch 22 calculation
            defaultParams.Catch22.BinSize = .1;%in seconds //100ms
            
            % Parameters for saving
            defaultParams.Save.Flag = false; %Flag if file should be save after analyses
            defaultParams.Save.Path = []; %Save path, if empty it uses the sorting path
            defaultParams.Save.Overwrite = false; %Flag if an existing file should be overwritten

        end

        function feat_table = getUnitFeatures(obj, unit_features)
            arguments
                obj
                unit_features string = ["ReferenceWaveform","ActivityFeatures"] %Alternatively WaveformFeatures
            end
            if unit_features == "all"
                unit_features = ["ActivityFeatures","WaveformFeatures","RegularityFeatures","Catch22","GraphFeatures"];
            end
            if class(obj) == "MEArecording"
                unit_array = [obj.Units];
            elseif class(obj) == "Unit"
                unit_array = obj;
            end
            feature_tables = {};
            for f = 1:length(unit_features)

                if unit_features(f) == "ReferenceWaveform"
                    aligned_wf = RecordingGroup.alignWaveforms(unit_array);
                    var_names = "Waveform" + (1:size(aligned_wf,1));
                    feature_tables{f} = array2table(aligned_wf',"VariableNames",var_names);

                elseif ismember(unit_features(f),["ACG","FullACG"])
                    acgs = [unit_array.(unit_features(f))];
                    norm_acgs = (acgs)./max(acgs);
                    norm_acgs(isnan(norm_acgs)) = 0;
                    var_names = unit_features(f) + (1:size(norm_acgs,1));
                    feature_tables{f} = array2table(norm_acgs',"VariableNames",var_names);
                else
                    feature_tables{f} = vertcat(unit_array.(unit_features(f)));
                end
            end
            feat_table = [feature_tables{:}];
        end
    end
    
    methods 

        function ParentPath = get.ParentRecordingPath(obj)
            path_parts = strsplit(obj.Metadata.InputPath,'/');
            ParentPath = strjoin([path_parts(1:end-1), "qc_output"],'/');
        end

        function acgs = get.ACGs(obj)
            ccgs = obj.Connectivity.CCG.CCGs;
            for u = 1:size(ccgs,1)
                acgs = nan(size(ccgs,[1,2]));
                for i = 1:size(ccgs,2)
                    acgs(:,i) = ccgs(:,i,i);
                end
            end
        end

        function acgs = get.FullACGs(obj)
            ccgs = obj.Connectivity.FullCCG.CCGs;
            for u = 1:size(ccgs,1)
                acgs = nan(size(ccgs,[1,2]));
                for i = 1:size(ccgs,2)
                    acgs(:,i) = ccgs(:,i,i);
                end
            end
        end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Analyses
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function unit_array = generateUnits(obj)
           [max_amplitudes, reference_electrode, norm_wf_matrix] = obj.generateWaveformMatrix();
           [firing_rates, unit_spike_times] = obj.calculateFiringRates(length(max_amplitudes));
           
           
           if isempty(obj.Parameters.QC.GoodUnits)
               good_units = performUnitQC(obj, max_amplitudes, firing_rates', unit_spike_times', norm_wf_matrix);
               obj.Parameters.QC.GoodUnits = good_units;
           else
               good_units = obj.Parameters.QC.GoodUnits;
               if length(good_units) ~= length(firing_rates)
                   error("Good unit IDs do not match. Did you use the correct ones?")
               end
               fprintf('Kept the %i good units provided\n', sum(good_units))
           end
           
           if sum(good_units) >= obj.Parameters.QC.N_Units
               good_amplitudes = max_amplitudes(good_units);
               good_unit_spike_times = unit_spike_times(good_units);
               good_wf_matrix = norm_wf_matrix(:,good_units);
               waveform_features = obj.inferWaveformFeatures(good_amplitudes, good_wf_matrix);
               good_template_ids = find(good_units);
               good_reference_electrodes = reference_electrode(good_units);
               unit_array = Unit();
               for u = 1:length(waveform_features)
                   unit_array(u) = Unit(obj, good_wf_matrix(:,u), good_reference_electrodes(u), good_unit_spike_times{u}, waveform_features{u});
                   unit_array(u).TemplateID = good_template_ids(u);
               end
               no_act_idx = arrayfun(@(x) isempty(x.ActivityFeatures),unit_array); 
               %If activity features were not computed due to not enough activity, we fill them to prevent issues when tracking units
%                act_table = vertcat(unit_array.ActivityFeatures);
%                empty_act_table = array2table(zeros(1,size(act_table,2)),'VariableNames',act_table.Properties.VariableNames);
                act_features = Unit.returnFeatureNames("act");
                empty_act_table = array2table(zeros(1,length(act_features)),'VariableNames',act_features);
               [unit_array(no_act_idx).ActivityFeatures] = deal(empty_act_table);
               if obj.Parameters.Analyses.Regularity
                   no_reg_idx = arrayfun(@(x) isempty(x.RegularityFeatures),unit_array);
                   %                    reg_table = vertcat(unit_array.RegularityFeatures);
                   %                    empty_reg_table = array2table(zeros(1,size(reg_table,2)),'VariableNames',reg_table.Properties.VariableNames);
                   reg_features = Unit.returnFeatureNames("reg");
                   empty_reg_table = array2table(zeros(1,length(reg_features)),'VariableNames',reg_features);
                   [unit_array(no_reg_idx).RegularityFeatures] = deal(empty_reg_table);
               end
               if obj.Parameters.Analyses.Catch22
                   no_c22_idx = arrayfun(@(x) isempty(x.Catch22),unit_array);
                   %                    c22_table = vertcat(unit_array.Catch22);
                   %                    empty_c22_table = array2table(zeros(1,size(c22_table,2)),'VariableNames',c22_table.Properties.VariableNames);
                   c22_features = Unit.returnFeatureNames("c22");
                   empty_c22_table = array2table(zeros(1,length(c22_features)),'VariableNames',"SC_" + c22_features);
                   [unit_array(no_c22_idx).Catch22] = deal(empty_c22_table);
               end
               
               obj.Units = unit_array;
               obj.updateSpikeTimes();
           else
               warning("Not enough good units found")
           end
       end
       
       function catch_22_table = run_catch_22(obj, spike_train, bin_size)
            arguments
                obj MEArecording
                spike_train (1,:) double = obj.Spikes.Times % in seconds
                bin_size (1,1) double = obj.Parameters.Catch22.BinSize
            end
            n_bins = round(obj.RecordingInfo.Duration / bin_size);
            binned = histcounts(spike_train,n_bins);
            [featureValues, featureNames] = catch22_all(binned');
            catch_22_table = array2table(featureValues','VariableNames',featureNames);
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
           Regularity.NetworkRegularityFrequency = F(freq_idx);
           Regularity.NetworkRegularityMagnitude = mag;
           
           norm_freq_domain = freq_domain/max(freq_domain);
           l = [];
           p = [];
           for i = 1:obj.Parameters.Regularity.N_peaks
               if length(norm_freq_domain) < 10
                   break
               else
                   [pks,locs,~,~] = findpeaks(norm_freq_domain,'NPeaks',1,'SortStr','descend');
                   norm_freq_domain = norm_freq_domain(locs:end);
                   l = [l; locs]; %#ok<AGROW>
                   p = [p; pks]; %#ok<AGROW>
                   if isempty(pks)
                       break
                   end
               end
           end
           log_p = log10(p)-min(log10(p)); % Make data positive to be able to fit
           try
               f = fit(cumsum(l),log_p,'exp1');
               Regularity.NetworkRegularityFit = f.b;
           catch
               Regularity.NetworkRegularityFit = NaN;
           end
           obj.NetworkFeatures.Regularity = struct2table(Regularity);
       end
       
       function inferOptimalBurstParameters(obj, plot_idx)
           arguments
               obj MEArecording
               plot_idx = false
           end
           [N, ISIn_ms, dunns_coeffs] = obj.inferISINparameters(); %"plot_idx",plot_idx
           ISI_N = ISIn_ms/1000;
           dunns_coeffs(ISI_N<0.01) = 0;
           [max_vals, max_idx] = maxk(dunns_coeffs,3);

           if ~isempty(max_vals)
               vals = sum(max_vals> max_vals(1)*0.9); %Check all Ns within 10% of the best N value
               candidate_idx = max_idx(1:vals);
           else
               candidate_idx = [];
           end
           % padded_N = [3 N N(end) + (N(end)-N(end-1))]; %Padding
           % candidate_Ns = [];%arrayfun(@(x) (padded_N(x-1)+1):(padded_N(x+1)-1),candidate_idx+1,'un',0); candidate_Ns = [candidate_Ns{:}];
           % if ~isempty(candidate_Ns)
           %     [refined_N, refined_ISIn_ms, refined_dunns_coeffs] = obj.inferISINparameters(N=candidate_Ns);
           %     [sorted,sort_idx] = sort(refined_dunns_coeffs,'descend');
           %     refined_vals = sum(sorted> sorted(1)*0.9);
           %     best_N = max(refined_N(sort_idx(1:refined_vals))); %If we have several good candidates, we take the higher N to be more conservative
           %     best_ISI_N = refined_ISIn_ms(refined_N == best_N)/1000;
           % else
           if isempty(candidate_idx)
               best_N = nan;
               best_ISI_N = nan;
           else
               best_N = max(N(candidate_idx));
               best_ISI_N = ISI_N(N==best_N);
           end
           obj.Bursts.best_N = best_N;
           obj.Bursts.N = N;
           obj.Bursts.best_ISI_N = best_ISI_N;
           obj.Bursts.ISI_N = ISI_N;
           obj.Bursts.dunns = dunns_coeffs;
       end

       function Burst = detectBursts(obj, ops)
           arguments
               obj MEArecording
               ops.N double = []
               ops.ISI_N double = []
               ops.merge_factor double = 0.5
               ops.plot_idx = false
           end
           if isempty(ops.N) | isempty(ops.ISI_N)
               if ~(isfield(obj.Bursts,'best_N') && isfield(obj.Bursts,'best_ISI_N'))
                   inferOptimalBurstParameters(obj, ops.plot_idx);
               end
               N = obj.Bursts.best_N;
               ISI_N = obj.Bursts.best_ISI_N;
           else
               N = ops.N;
               ISI_N = ops.ISI_N;
           end

           if isnan(obj.Bursts.best_N) %If no good parameters could be inferred
               obj.Bursts.T_start = [];
               obj.Bursts.T_end = [];
           else
               Burst = detectBurstsISIN(obj, N, ISI_N);
               obj.Bursts.T_start = Burst.T_start;
               obj.Bursts.T_end = Burst.T_end;
               obj.pruneBursts();
               obj.mergeBursts(ops.merge_factor);
           end
           
           if ops.plot_idx
               obj.PlotBurstCheck(); title(sprintf("N=%i; ISI_N=%.3f, Units=%i",N,ISI_N,length(obj.Units)))
           end
       end

       function mergeBursts(obj,merge_factor, input)
           arguments
               obj MEArecording
               merge_factor (1,1) double = 0.5 %Bursts with an IBI of merge_factor * best_ISI_N will be merged
               input (1,1) string = "processed" %"raw"
           end
           if input== "raw" && isfield(obj.Bursts,'raw')
               Burst = obj.Bursts.raw;
           else
               Burst = obj.Bursts;
               obj.Bursts.raw = obj.Bursts;
           end
           
           N_bursts = length(Burst.T_start);
           IBI_merge = merge_factor * obj.Bursts.best_ISI_N;
           if N_bursts > 1
               %            Merge bursts with too short IBIs (< IBI_merge)
               short_burst_idx = find(Burst.T_start(2:end) - Burst.T_end(1:end-1)> IBI_merge);
               obj.Bursts.T_start = Burst.T_start([1 short_burst_idx+1]);
               obj.Bursts.T_end = Burst.T_end([short_burst_idx length(Burst.T_end)]);
               fprintf('Merged %i bursts\n', N_bursts - length(obj.Bursts.T_start))
           end

       end

       function [clean_spikes, clean_units] = remove_tonic_units(obj,ops)
           arguments
               obj MEArecording
               ops.plot_s = 10
               ops.plot_idx = false
           end
           act_table = obj.getUnitFeatures(obj.Units,"ActivityFeatures");
           tonic_ids = find((isoutlier(act_table.FiringRate,"ThresholdFactor",1) & act_table.FiringRate > median(act_table.FiringRate)) & ...
               act_table.CVInterSpikeInterval < mean(act_table.CVInterSpikeInterval));
           spike_times = obj.Spikes.Times;
           spike_units = obj.Spikes.Units;
           clean_idx = ~ismember(spike_units,tonic_ids);
           clean_units = findgroups(spike_units(clean_idx));
           clean_spikes = spike_times(clean_idx);
           tonic_spikes = spike_times(~clean_idx);
           tonic_units = findgroups(spike_units(~clean_idx));
           if ops.plot_idx
               figure('Color','W');
               subplot(121); scatter(clean_spikes(clean_spikes < ops.plot_s ),clean_units(clean_spikes < ops.plot_s ),0.1,'k.');title('Clean')
               subplot(122); scatter(tonic_spikes(tonic_spikes < ops.plot_s ), tonic_units(tonic_spikes < ops.plot_s ),0.1,'k.');title('Tonic')
           end
       end

       function [N_array, min_isis, dunns_coeffs] = inferISINparameters(obj, ops)
           arguments
               obj MEArecording
               ops.N (1,:) double = []
               ops.Steps = 10.^(-3:.025:2)
               ops.plot_idx = false
           end
           [min_isis, N_array, prom] = obj.calculateMinISIN(N=ops.N, Steps=ops.Steps, plot_idx=ops.plot_idx);
            N_array = N_array(min_isis > 0); % Remove Ns for which we did not see a second peak in the histogram
            min_isis(min_isis == 0) = [];
           [SpikeTimes,~] = obj.remove_tonic_units(plot_idx=ops.plot_idx); 
           
           dunns_coeffs = zeros(size(N_array));
           N_spikes = 10000; %Cutoff to prevent the distance matrix from becoming too big
           for ni = 1:length(N_array)
               i = 0; max_i = max([1, floor(length(SpikeTimes)/N_spikes)]); ISI_N = -1;

               while min(ISI_N) > min_isis(ni) || max(ISI_N) < min_isis(ni) && i<max_i %iterate through the spike times to find suitable cutout
                   spike_cutout = SpikeTimes((i*N_spikes)+1:min([(i+1)*N_spikes,length(SpikeTimes)]));
                   ISI_N = 1e3 * (spike_cutout(N_array(ni):end) - spike_cutout(1:end-(N_array(ni)-1)));
                   i = i+1;
               end

               if min(ISI_N) > min_isis(ni) || max(ISI_N) < min_isis(ni)
                   dunns_coeffs(ni) = 0;
               else
                   idx = ones(size(ISI_N));
                   idx(ISI_N>min_isis(ni)) = 2;
                   distM=squareform(pdist(ISI_N));
                   dunns_coeffs(ni) = dunns(2, distM, idx);
               end
           end
       end

       function [min_isis, N_array, prom] = calculateMinISIN(obj, ops)
           arguments
               obj MEArecording
               ops.N (1,:) double = []
               ops.Steps = 10.^(-3:.025:2)
               ops.plot_idx = false
           end
           [SpikeTimes, SpikeUnits] = obj.remove_tonic_units();

           if isempty(ops.N)
               min_N = max([4 ceil(length(unique(SpikeUnits))/30)]);
               max_N = ceil(length(unique(SpikeUnits))/2);
               if (max_N - min_N) < 20
                   N_array = min_N:max_N;
               else
                   N_array = round(linspace(min_N,max_N,20));
                   % N_array = min_N:max_N;
               end
           else
               N_array = ops.N;
           end
           
           prom = zeros(1,length(N_array));
           min_isis = zeros(1,length(N_array));
           if ops.plot_idx
               figure;tiledlayout('flow');hold on;colormap(hsv)
           end
           for i = 1:length(N_array)
               FRnum = max([ 2 N_array(i)]);
               ISI_N = SpikeTimes( FRnum:end ) - SpikeTimes( 1:end-(FRnum-1) );
               n = histc( ISI_N*1000, ops.Steps*1000 );
               n = smooth( n, 'lowess' );
               [pks,b,w,p] = findpeaks(n);
               if length(b) < 2
                   min_isis(i) = 0;
               else
                   [sorted_p,sorted_idx] = sort(p,'descend');
                   b = b(sorted_idx(1:2));
                   [~,m_idx] = min(n(min(b):max(b)));
                   isi_idx = (m_idx+min(b)-1);
                   min_isis(i) = ops.Steps(isi_idx)*1000;
                   prom(i) = sorted_p(2);
               end
               if ops.plot_idx
                   nexttile
                   plot(ops.Steps*1000,n/sum(n),'--');xscale('log');title(sprintf('N %i, ISI_N %.3f',N_array(i),min_isis(i)))
               end
           end
       end

       function pruneBursts(obj,ops)
           arguments
               obj MEArecording
               ops.bins = 200
               ops.min_width = 0.005
           end

           pruned.T_start = nan(size(obj.Bursts.T_start));
           pruned.T_end = nan(size(obj.Bursts.T_end));
           % figure('Color','W');tiledlayout('flow','TileSpacing','compact','Padding','compact')
           for iburst = 1:length(obj.Bursts.T_start)
               padding = max([obj.Bursts.best_ISI_N, (obj.Bursts.T_end(iburst) - obj.Bursts.T_start(iburst)) * 0.5]);
               padded_start = obj.Bursts.T_start(iburst) - padding;
               padded_end = obj.Bursts.T_end(iburst) + padding;
               if ops.bins > 1 % interpret as num_bins if > 1
                   bin_width = max([ops.min_width, (padded_end - padded_start) / ops.bins]);
               else
                   bin_width = ops.bins;
               end
               burst_activity = obj.Spikes.Times(obj.Spikes.Times <= padded_end & obj.Spikes.Times >= padded_start);
               binned_activity = histcounts(burst_activity,padded_start:bin_width:padded_end);
               smoothed_activity = smoothdata(binned_activity,'lowess','SmoothingFactor',0.5) + 1;
               cum_activity_l = cumsum(smoothed_activity);
               cum_activity_r = cumsum(smoothed_activity(end:-1:1));
               threshold_l = triangle_threshold(cum_activity_l, 'L', false);
               threshold_r = triangle_threshold(cum_activity_r(end:-1:1), 'R', false);

               % Sanity checks
               checks(1) = mean(binned_activity(threshold_l:threshold_r)) > mean(binned_activity);
               checks(2) = threshold_l > 1 && threshold_r < length(binned_activity);
               checks(3) = threshold_l < threshold_r;
               % nexttile
               % plot(binned_activity);xline([threshold_l,threshold_r],'r')
               % axis off 
% pruned_activity = obj.Spikes.Times(obj.Spikes.Times > (padded_start + threshold_l * bin_width) & obj.Spikes.Times < padded_end - (length(smoothed_activity) - threshold_r) * bin_width);
% pruned_binned = histcounts(pruned_activity,(padded_start + threshold_l * bin_width):bin_width:(padded_end - (length(smoothed_activity) - threshold_r) * bin_width));
               if all(checks)
                   pruned.T_start(iburst) = padded_start + threshold_l * bin_width;
                   pruned.T_end(iburst) = padded_end - (length(smoothed_activity) - threshold_r) * bin_width;
               else
                   pruned.T_start(iburst) = obj.Bursts.T_start(iburst);
                   pruned.T_end(iburst) = obj.Bursts.T_end(iburst);
                   % title('failed')
               end
           end
           obj.Bursts.T_start = pruned.T_start;
           obj.Bursts.T_end = pruned.T_end;
       end

       function [inburst_activity, inburst_units] =  getBurstStatistics(obj, ops)
           arguments
               obj MEArecording
               ops.binning = 0.01
           end
           if isempty(obj.Bursts) %Burst detection was not performed
               detectBursts(obj, merge_factor=obj.Parameters.Bursts.MergeFactor);
           end
           if length(obj.Bursts.T_start) < 3 %Not enough bursts were detected
               Burst.MeanInterBurstInterval = 0;
               Burst.BurstDuration = 0;
               Burst.RiseTime = 0;
               Burst.FallTime = 0;
               Burst.PeakFR = 0;
               Burst.StdFR = 0;
               Burst.StdBurstDuration = 0;
               Burst.StdInterBurstInterval = 0;
               Burst.IntraFiringRate = 0;
               Burst.InterFiringRate = 0;
           else
               IBIs = obj.Bursts.T_start(2:end) - obj.Bursts.T_end(1:end-1);
               BDs = obj.Bursts.T_end - obj.Bursts.T_start;
               
               [RiseTime, FallTime, PeakFR, StdFR] = deal(nan(1,length(BDs)));
               Fs = round(1/ops.binning);
               inburst_activity = cell(1,length(BDs));
               inburst_units = cell(1,length(BDs));
               for iburst = 1:length(BDs)
                   try
                       burst_activity = obj.Spikes.Times(obj.Spikes.Times <= obj.Bursts.T_end(iburst) & obj.Spikes.Times >= obj.Bursts.T_start(iburst));
                       binned_activity = histcounts(burst_activity,'BinWidth',ops.binning);
                       [max_binned, imax_binned] = max(binned_activity);
                       norm_activity = smoothdata(binned_activity/max_binned,'lowess');
                       RiseTime(iburst) = imax_binned * ops.binning;
                       FallTime(iburst) = (length(binned_activity) - imax_binned) * ops.binning;
                       %%% Old approch, but seems too unreliable %%%

                       % if length(norm_activity) < 4
                       %     [peak_1,peak_1_idx] = max(norm_activity);
                       %     peak_2 = peak_1; peak_2_idx = peak_1_idx;
                       % else
                       %     [peak_1,peak_1_idx] = findpeaks(norm_activity,"MinPeakProminence",0.3,NPeaks=1);
                       %     [peak_2,peak_2_idx] = findpeaks(norm_activity(end:-1:1),"MinPeakProminence",0.3,NPeaks=1);
                       %     peak_2_idx = length(norm_activity) - peak_2_idx + 1;
                       % end
                       % if isempty(peak_1)
                       %     [peak_1,peak_1_idx] = max(norm_activity);
                       % end
                       % if isempty(peak_2)
                       %     [peak_2,peak_2_idx] = max(norm_activity);
                       % end
                       % 
                       % pre_lower_level = min(norm_activity(1:peak_1_idx));
                       % post_lower_level = min(norm_activity(peak_2_idx:end));
                       % % Using mean to account for the rare case where rise/fall is split
                       % RiseTime(iburst) = mean(risetime([ones(1,Fs)*pre_lower_level norm_activity(1:peak_1_idx) ones(1,Fs)*peak_1],Fs));
                       % FallTime(iburst) = mean(falltime([ones(1,Fs)*peak_2 norm_activity(peak_2_idx:end) ones(1,Fs)*post_lower_level],Fs));

                       %%% Old approch, but seems too unreliable %%%

                       PeakFR(iburst) = max_binned./length(obj.Units) * Fs;
                       StdFR(iburst) = std(norm_activity);
                   catch ME
                       warning(ME.message)
                   end
                   
               end
               [inburst, outburst] = obj.getBurstSpikes();
               inter_spikes = cellfun(@length, outburst.spikes);
               inter_fr = inter_spikes(2:end-1) ./ (IBIs * length(obj.Units));
               intra_spikes = cellfun(@length, inburst.spikes);
               intra_fr = intra_spikes ./ (BDs * length(obj.Units));

               cellfun_input = {IBIs, BDs, RiseTime, FallTime, PeakFR, StdFR};
               
               % Perform outlier removal if method is specified
               if ~isempty(obj.Parameters.Outlier.Method)
                   cellfun_input = cellfun(@(x) rmoutliers(x,obj.Parameters.Outlier.Method,'ThresholdFactor',obj.Parameters.Outlier.ThresholdFactor),...
                       cellfun_input,'un',0);
               end
               feature_means = cellfun(@nanmedian, cellfun_input,'un',0);
               
               [Burst.MeanInterBurstInterval,...
                   Burst.BurstDuration,...
                   Burst.RiseTime,...
                   Burst.FallTime,...
                   Burst.PeakFR,...
                   Burst.StdFR] = feature_means{:};
               Burst.StdBurstDuration = std(BDs);
               Burst.StdInterBurstInterval = std(IBIs);
               Burst.IntraFiringRate = median(intra_fr);
               Burst.InterFiringRate = median(inter_fr);
           end
           obj.NetworkFeatures.Burst = struct2table(Burst);
       end

       function [inburst, outburst] = getBurstSpikes(obj)
           inburst.spikes = cell(1,length(obj.Bursts.T_start));
           inburst.units = cell(1,length(obj.Bursts.T_start));
           outburst.spikes = cell(1,length(obj.Bursts.T_start)+1);
           outburst.units = cell(1,length(obj.Bursts.T_start)+1);

           for iburst = 1:length(obj.Bursts.T_start)
               burst_idx = obj.Spikes.Times > obj.Bursts.T_start(iburst) & obj.Spikes.Times < obj.Bursts.T_end(iburst);
               inburst.spikes{iburst} = obj.Spikes.Times(burst_idx);
               inburst.units{iburst} = obj.Spikes.Units(burst_idx);
               if iburst == 1
                   outburst_idx = obj.Spikes.Times > 0 & obj.Spikes.Times < obj.Bursts.T_start(iburst);
               else
                   outburst_idx = obj.Spikes.Times > obj.Bursts.T_end(iburst-1) & obj.Spikes.Times < obj.Bursts.T_start(iburst);
               end
               outburst.spikes{iburst} = obj.Spikes.Times(outburst_idx);
               outburst.units{iburst} = obj.Spikes.Units(outburst_idx);
           end
           outburst_idx = obj.Spikes.Times > obj.Bursts.T_end(iburst);
           outburst.spikes{iburst+1} = obj.Spikes.Times(outburst_idx);
           outburst.units{iburst+1} = obj.Spikes.Units(outburst_idx);
       end

       function inferPowerSpectrum(obj, ops)
           arguments
               obj MEArecording
               ops.binning double = 0.001 %in seconds
               ops.min_length double = 0 %in seconds
               ops.min_participation double = 0.2 %minimum ratio of bursts a unit is in to be considered for analysis
               ops.concatenated logical = true %Flag whether bursts should be concatenated or not
               ops.overwrite logical = false
               % only relevant if concatenated=true
               ops.duration double = 30 %we split up the concatenated signal into chunks of this length
               ops.overlap double = 0.1 %with this much overlap
               ops.freq_range (1,2) double = [1, 500]
           end

            if isfield(obj.Bursts,'FullPowerSpectrum') && ~ops.overwrite
                disp('Power spectrum already computed.')
                return
            end

           [inburst, outburst] = obj.getBurstSpikes();
           inburst_data = generate_fooof_input(inburst.spikes,{},binning=ops.binning,min_length=ops.min_length,...
               min_participation=ops.min_participation,concatenated=ops.concatenated,duration=ops.duration);
           inburst_freq = calc_pow(inburst_data,freq_range=ops.freq_range);
           outburst_data = generate_fooof_input(outburst.spikes,{},binning=ops.binning,min_length=ops.min_length,...
               min_participation=ops.min_participation,concatenated=true,duration=ops.duration);
           outburst_freq = calc_pow(outburst_data,freq_range=ops.freq_range);

           obj.Bursts.FullPowerSpectrum.Inburst = inburst_freq;
           obj.Bursts.FullPowerSpectrum.Outburst = outburst_freq;
       end

       function inferFooofSpectrum(obj, ops)
           arguments
               obj MEArecording
               ops.binning double = 0.001 %in seconds
               ops.min_length double = 0 %in seconds
               ops.min_participation double = 0.2 %minimum ratio of bursts a unit is in to be considered for analysis
               ops.concatenated logical = true %Flag whether bursts should be concatenated or not
               ops.overwrite logical = false
               % only relevant if concatenated=true
               ops.duration double = 30 %we split up the concatenated signal into chunks of this length
               ops.overlap double = 0.1 %with this much overlap
               ops.freq_range (1,2) double = [10, 70]
           end
           if isfield(obj.Bursts,'PowerSpectrum') && ~ops.overwrite
                disp('Power spectrum already computed.')
                return
           end
           [inburst, outburst] = obj.getBurstSpikes();

           inburst_data = generate_fooof_input(inburst.spikes,{},binning=ops.binning,min_length=ops.min_length,...
               min_participation=ops.min_participation,concatenated=ops.concatenated,duration=ops.duration);

           [fooof_inburst.oscillatory,fooof_inburst.fractal, fooof_inburst.original] = foof(inburst_data,freq_range = ops.freq_range);
           obj.Bursts.PowerSpectrum.Inburst = fooof_inburst;

           outburst_data = generate_fooof_input(outburst.spikes,{},binning=ops.binning,min_length=ops.min_length,...
               min_participation=ops.min_participation,concatenated=true,duration=ops.duration);
           [fooof_outburst.oscillatory,fooof_outburst.fractal, fooof_outburst.original] = foof(outburst_data,freq_range = ops.freq_range);
           obj.Bursts.PowerSpectrum.Outburst = fooof_outburst;
       end

       function PlotFOOOFresults(obj)
           if ~isfield(obj.Bursts,'PowerSpectrum')
               obj.inferPowerSpetrum();
           end
           ib = obj.Bursts.PowerSpectrum.Inburst;
           ob = obj.Bursts.PowerSpectrum.Outburst;

           figure('Color','W'); tiledlayout(1,3)
           nexttile
           plot(ib.oscillatory.freq,log(ib.oscillatory.powspctrm))
           hold on
           plot(ob.oscillatory.freq,log(ob.oscillatory.powspctrm))

           nexttile
           plot(log(ib.fractal.freq),log(ib.fractal.powspctrm))
           hold on
           plot(log(ob.fractal.freq),log(ob.fractal.powspctrm))

           nexttile
           plot(log(ib.original.freq),log(ib.original.powspctrm))
           hold on
           plot(log(ob.original.freq),log(ob.original.powspctrm))
           legend(["Inburst","Outburst"])
           sgtitle(obj.Metadata.ChipID + " " + obj.Metadata.AAV + " " + obj.Metadata.Concentration)


       end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Connectivity functions
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function inferConnectivity(obj, alg)
          arguments
             obj MEArecording
             alg (1,:) string %CCG (English et al.), DDC ( et al.), STTC (Cutts et al.) // can also be several 
          end
          
          for a = 1:length(alg)
              if alg(a) == "FullCCG"
                  obj.inferConnectivityCCG("FullCCG");
              else
                  fh = str2func("inferConnectivity" + alg(a));
                  fh(obj);
              end
          end
       end

       function [ccgs,t] = calculateFullCCG(obj,args) %Use only if the recording was split from a longer/concatenated recording
           arguments
               obj MEArecording
               args.binsize (1,1) double = 0.0005
               args.duration (1,1) double = 0.1
           end
           spike_times = readNPY(fullfile(obj.ParentRecordingPath,'spike_times.npy'));
           spike_templates = readNPY(fullfile(obj.ParentRecordingPath,'spike_templates.npy'));
           good_ids = [obj.Units.TemplateID];
           spike_times = double(spike_times(ismember(spike_templates,good_ids)))/obj.RecordingInfo.SamplingRate;
           spike_templates = spike_templates(ismember(spike_templates,good_ids));
           spike_templates = findgroups(spike_templates);

           [ccgs,t] = CCG(spike_times,spike_templates,'Fs',(1/obj.RecordingInfo.SamplingRate),'binsize',args.binsize,'duration',args.duration);
           obj.Connectivity.FullCCG.CCGs = ccgs;
           obj.Connectivity.FullCCG.args = args;
       end

       function [ccg,t] = calculateCCG(obj)
           if isempty(which('CCGHeart'))
               cd(fullfile(obj.getParentPath(),'Functions'))
               mex('CCGHeart.c')
           end
           [ccg,t] = CCG(obj.Spikes.Times,obj.Spikes.Units,'binSize',obj.Parameters.CCG.BinSize,'duration',obj.Parameters.CCG.Duration,'Fs',1/obj.RecordingInfo.SamplingRate);
           disp('Finished CCG calculations')
       end

       function ccgs = assembleCCGs(obj, rg, culture_idx)
           arguments
               obj MEArecording
               rg RecordingGroup
               culture_idx (1,:) double = []
           end
           culture = rg.findCulture(obj);
           if ~isempty(culture_idx)
               culture = culture(culture_idx);
           end
           ccgs = zeros(size(culture(1).Connectivity.CCG.CCGs));
           for r = 1:length(culture)
               ccgs = ccgs + culture(r).Connectivity.CCG.CCGs;
           end
       end

       function results = inferConnectivityGLM(obj)
           addpath(genpath('/home/phornauer/Git/extended-GLM-for-synapse-detection'))
           hyperparameter = obj.Parameters.GLM;
           try
               ccgs = obj.Connectivity.FullCCG.CCGs;
           catch
               [ccgs,t] = calculateCCG(obj);
               obj.Connectivity.CCG.CCGs = ccgs;
           end
           ref_els = [obj.Units.ReferenceElectrode];
           location.x = obj.RecordingInfo.ElectrodeCoordinates(ref_els,1);
           location.y = obj.RecordingInfo.ElectrodeCoordinates(ref_els,2);
           NN = size(ccgs,2);
           CCG = cell(NN);

           distance = zeros(NN,NN);
           ignore = zeros(NN,NN);
           for i = 1:NN
               for j = 1:NN
                   distance(i,j) = sqrt((location.x(i)-location.x(j))^2 + (location.y(i)-location.y(j))^2);
                   ccg = ccgs(:,i,j);
                   CCG{i,j} = ccg';
                   if sum(ccgs(:,i,j)) < hyperparameter.min_spikes || i==j
                       ignore(i,j) = 1;
                   end
               end
           end
           X = learning_basis(CCG,ignore);
           if isinparfor
               for pre = 1:NN % use parfor for parallel computing
                   model_fits(pre) = extendedGLM(CCG(pre,:), X, distance(pre,:),hyperparameter,ignore(pre,:));
               end
           else
               parfor pre = 1:NN % use parfor for parallel computing
                   model_fits(pre) = extendedGLM(CCG(pre,:), X, distance(pre,:),hyperparameter,ignore(pre,:));
               end
           end
           results = detect_cnx(model_fits,ignore,hyperparameter.threshold);
           obj.Connectivity.concatenated.GLM = results;
       end

       function inferConnectivityCCG(obj,ccg_type)
           arguments
            obj MEArecording
            ccg_type string = "CCG" %or "FullCCG"
           end
           if ~isfield(obj.Connectivity,ccg_type) || ~isfield(obj.Connectivity.(ccg_type),'CCGs')
               fh = str2func("calculate" + ccg_type);
               fh(obj);
           end
           [ccgR1,tR] = obj.calculateCCG();
           if size(ccgR1,2) < length(obj.Units) %Pad if the last unit(s) have no spikes
               padded_ccg = zeros(size(ccgR1,1),length(obj.Units),length(obj.Units));
               padded_ccg(:,1:size(ccgR1,2),1:size(ccgR1,2)) = ccgR1;
               ccgR1 = padded_ccg;
           end
           binSize = obj.Parameters.CCG.BinSize; %.5ms
           conv_w = obj.Parameters.CCG.Conv_w;  % 10ms window
           alpha = obj.Parameters.CCG.Alpha; %high frequency cut off, must be .001 for causal p-value matrix
%            Fs = 1/obj.RecordingInfo.SamplingRate;
           
           nCel=size(ccgR1,2);
           
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
           
%            silent_units = find(arrayfun(@(x) isempty(x.SpikeTimes),obj.Units));
           
           obj.Connectivity.(ccg_type).ExcitatoryConnection = sig_con;
           obj.Connectivity.(ccg_type).InhibitoryConnection = sig_con_inh;
           obj.Connectivity.(ccg_type).bd = con_mat;
           obj.Connectivity.(ccg_type).CCGs = ccgR;
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
%            V_obs = smoothdata(V_obs,'gaussian',3);
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
           obj.Connectivity.DDC.wu = dCov * inv(B);
       end
       
       function inferConnectivitySTTC(obj)
           if ~isempty(obj.Spikes.Times)
               empirical_sttc = obj.empiricalSTTCmatrix();
               surrogate_mat= nan([size(empirical_sttc), obj.Parameters.STTC.N_Surrogates]);
               for i = 1:obj.Parameters.STTC.N_Surrogates
                   surrogate_mat(:,:,i) = obj.surrogateSTTCmatrix();
               end
               STTC.threshold = prctile(surrogate_mat,obj.Parameters.STTC.Percentile,3);
               STTC.wu = empirical_sttc;
               STTC.bu = STTC.wu > STTC.threshold;
               obj.Connectivity.STTC = STTC;
           else
               STTC.wu = zeros(length(obj.Units));
               STTC.bu = STTC.wu;
               obj.Connectivity.STTC = STTC;
           end
       end
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Graph features
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function inferGraphFeatures(obj, alg)
          arguments
             obj MEArecording
             alg (1,:) string = [] %Connectivity inference algorithms to compute graphs for
          end
          
          if isempty(alg) %Compute features for all available inferred graphs
              alg = string(fields(obj.Connectivity));
              assert(~isempty(alg),"No connectivity inference results found")
          end
          
          full_nw_table = table();
          full_unit_table = table();
          for a = 1:length(alg)
             assert(isfield(obj.Connectivity,alg(a)),"No connectivity inference results found")
             if isfield(obj.Connectivity.(alg(a)),"bd") %Binary directed graphs
                 con_mat = obj.Connectivity.(alg(a)).bd;
                 density = density_dir(con_mat);
                 clustering_coef = clustering_coef_bd(con_mat);
                 assortativity = assortativity_bin(con_mat,1);
%                  rich_club = rich_club_bd(con_mat,5);
                 nw_feat = ["Density","Assortativity"];%,"RichClub" + (1:5)];
                 nw_val = [density, assortativity];%, rich_club];
                 
             elseif isfield(obj.Connectivity.(alg(a)),"bu") %Binary undirected graphs
                 con_mat = obj.Connectivity.(alg(a)).bu;
                 density = density_und(con_mat);
%                  XY = obj.RecordingInfo.ElectrodeCoordinates([obj.Units.ReferenceElectrode],:);
%                  n = 1000;
%                  tol = 1e-6;
%                  [N,E] = rentian_scaling_2d(con_mat,XY,n,tol);
%                  rm_idx = N == 0 | E == 0;
%                  N(rm_idx) = []; E(rm_idx) = [];
%                  try
%                      [b,~] = robustfit(log10(N),log10(E));
%                      rents_exponent = b(2,1);
%                  catch
rents_exponent = 0; 
%                  end
                 clustering_coef = clustering_coef_bu(con_mat);
                 assortativity = assortativity_bin(con_mat,0);
%                  rich_club = rich_club_bu(con_mat,5);
                 nw_feat = ["Density","RentExponent","Assortativity"];%,"RichClub" + (1:5)];
                 nw_val = [density, rents_exponent, assortativity];%, rich_club];
                 
             else 
                 disp("No binary connectivity matrix available")
             end
             
             global_efficiency = efficiency_bin(con_mat);
             local_efficiency = efficiency_bin(con_mat,1);
%              [s,S] = motif3struct_bin(con_mat);
%              s = s./nchoosek(length(con_mat),3); %Normalize
%              S = S./(nchoosek(length(con_mat),3) - nchoosek(length(con_mat) - 1,3)); %Normalize
%              [f,F] = motif3funct_bin(con_mat);
%              f = f./nchoosek(length(con_mat),3); %Normalize
%              F = F./(nchoosek(length(con_mat),3) - nchoosek(length(con_mat) - 1,3)); %Normalize
             eigen_centrality = eigenvector_centrality_und(double(con_mat));
             betweenness = betweenness_bin(con_mat);
             nw_feat = [nw_feat, "GlobalEfficiency"] + "_" + alg(a);%, "StructuralMotif" + (1:13), "FunctionalMotif" + (1:13)] + "_" + alg(a);
             nw_val = [nw_val, global_efficiency];%, s', f'];
             nw_val(isnan(nw_val) | isinf(nw_val)) = 0;
             nw_table = array2table(nw_val,"VariableNames",nw_feat);
             unit_feat = ["ClusteringCoefficient","LocalEfficiency","EigenCentrality","Betweenness"] + "_" + alg(a);%,"UnitStructuralMotif" + (1:13), "UnitFunctionalMotif" + (1:13)];
             unit_val = [clustering_coef, local_efficiency, eigen_centrality, betweenness'];%, S', F'];
             unit_val(isnan(unit_val) | isinf(unit_val)) = 0;
             unit_table = array2table(unit_val,"VariableNames",unit_feat);
             full_nw_table = [full_nw_table nw_table];
             full_unit_table = [full_unit_table unit_table];
          end
          obj.NetworkFeatures.GraphFeatures = full_nw_table;
          for u = 1:size(full_unit_table,1)
              obj.Units(u).GraphFeatures = full_unit_table(u,:);
          end
       end
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Feature value retrieval
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function merged_table = concatenateClusteredFeatures(obj,cluster_id, feature_group)
           arguments
               obj MEArecording
               cluster_id (1,:) {isnumeric} = 0 %Cluster IDs to be concatenated; 0 concatenates all
               feature_group string = "all" %"ActivityFeatures","WaveformFeatures" or "all"
           end
           for rec = 1:length(obj)
               unit_array = [obj(rec).Units];
               if ~isempty([unit_array.ClusterID])
                   if cluster_id == 0
                       cluster_id = 1:obj(rec).NumUnitClusters;
                   end
                   
                   feat_array = [obj(rec).ClusteredFeatures];
                   fnames = fieldnames([feat_array{:}]);
                   
                   clust_table = cell(1,max(cluster_id));
                   clust_cell = struct2cell([feat_array{:}]);
                   for i = cluster_id
                       if feature_group == "all"
                           feat_idx = 1:length(fnames);
                       else
                           feat_idx = find(contains(fnames, feature_group));
                       end
                       
                       clust_table{i} = [clust_cell{feat_idx,:,i}];
                       clust_table{i}.Properties.VariableNames = string(clust_table{i}.Properties.VariableNames) + "_" + num2str(i);
                   end
                   full_table = [clust_table{:}];
               else
                   error("Unit clustering needs to be performed before using clustered features")
               end
               merged_table(rec,:) = full_table;
           end
       end
       
       function feat_table = getRecordingFeatures(obj,network_features, unit_features, useClustered) %Usable for unit and network features
           arguments
               obj MEArecording
               network_features string = "all" %"all" or []
               unit_features string = "all" %"ActivityFeatures","WaveformFeatures" or "all"
               useClustered logical = false
           end
           if ~isempty(network_features)
               fname_array = arrayfun(@(x) fieldnames(x.NetworkFeatures), obj,'un',0);
               [gc, fnames] = groupcounts(vertcat(fname_array{:}));%groupcounts shuffles fieldnames so we have to extract the original order again
               rm_field = string(fnames{gc ~= length(fname_array)});
               fnames = fnames(gc == length(fname_array));
               % fnames = fieldnames([obj.NetworkFeatures]);
               for i = 1:length(obj)
                   network_struct = obj(i).NetworkFeatures;
                   if isfield(network_struct, rm_field)
                    network_struct = rmfield(network_struct, rm_field);
                   end
                   network_array(i) = network_struct;
               end
                fnames = fieldnames(network_array);
               if network_features == "all"
                   nw_idx = 1:length(fnames);
               else
                   nw_idx = find(contains(fnames,network_features));
               end
               
               network_cell = squeeze(struct2cell(network_array));
               network_cell = reshape(network_cell,[],length(obj)); %Ensure correct orientation of cell array
               for i = 1:size(network_cell,2)
                  network_table(i,:) = [network_cell{nw_idx,i}]; 
               end
           else
               network_table = table();
           end
           
           if ~isempty(unit_features)
               if useClustered
                   unit_table = concatenateClusteredFeatures(obj, 0, unit_features);
               else
                   % fname_array = arrayfun(@(x) fieldnames(x.UnitFeatures), obj,'un',0);
                   % [gc, fnames] = groupcounts(vertcat(fname_array{:}));
                   % fnames = fnames(gc == length(fname_array));
                   fnames = fieldnames([obj.UnitFeatures]);
                   
                   if unit_features == "all"
                       feat_idx = 1:length(fnames);
                   else
                       feat_idx = find(contains(fnames, unit_features));
                       if length(feat_idx) ~= length(unit_features)
                           error(join('Could not find: ' + unit_features(~contains(unit_features,fnames))))
                       end
                   end
                   unit_cell = squeeze(struct2cell([obj.UnitFeatures]));
                   unit_cell = reshape(unit_cell,[],length(obj)); %Ensure correct orientation of cell array
                   for i = 1:size(unit_cell,2)
                       unit_table(i,:) = [unit_cell{feat_idx,i}];
                   end
               end
           else
               unit_table = table();
           end
           
           feat_table = [network_table unit_table];
       end
       
      
       function close_units = findCloseUnits(obj, unit_id, min_dist, max_dist)
           arguments
               obj MEArecording
               unit_id
               min_dist = 50
               max_dist = 150
           end
           coor = obj.getElectrodeCoordinates();
           ref_els = [obj.Units.ReferenceElectrode];
           ref_coor = coor(ref_els,:);
           distances = squareform(pdist(ref_coor));
           unit_dist = distances(unit_id,:);
           close_units = find(unit_dist < max_dist & unit_dist > min_dist);
           close_units(close_units==unit_id) = [];
       end
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Plots
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function PlotNetworkScatter(obj, time_cutout, color)
           arguments
               obj MEArecording
               time_cutout (1,2) double = [0 180]
               color = 'k'
           end
           %           figure('Color','w');
           spikes = obj.Spikes.Times(obj.Spikes.Times > time_cutout(1) & obj.Spikes.Times < time_cutout(2));
           units = obj.Spikes.Units(obj.Spikes.Times > time_cutout(1) & obj.Spikes.Times < time_cutout(2));
           plot(spikes, units, '.','Color',color,'MarkerSize', 0.1)
           xlabel('Time [s]')
           ylabel('Unit ID')
           xticks(linspace(time_cutout(1),time_cutout(2),5))
%            xticklabels(unique([xt time_cutout]))
           axis tight
           box off
           xlim(time_cutout)
       end
       
       function PlotNetworkScatterHistogram(obj, time_cutout, binning, color)
           arguments
               obj MEArecording
               time_cutout (1,2) double = [0 180],
               binning = obj.Parameters.Bursts.Binning
               color = 'k'
           end
           spikes = obj.Spikes.Times(obj.Spikes.Times > time_cutout(1) & obj.Spikes.Times < time_cutout(2));
           units = obj.Spikes.Units(obj.Spikes.Times > time_cutout(1) & obj.Spikes.Times < time_cutout(2));
           duration = [time_cutout(2) - time_cutout(1)];
           N_bins = [round(duration / binning); length(obj.Units)];
%            figure('Color','w');
           scatterhistogram(spikes,units,'MarkerSize',1,'Color',color,'HistogramDisplayStyle','stairs','NumBins', N_bins);
           
       end

       function PlotBurstCheck(obj,type,ops)
           arguments
               obj MEArecording
               type string = "hist" %"scatter"
               ops.cutout = [0,50]
               ops.binning = 0.01
               ops.input = "processed" %"raw" to use unmerged+untrimmed data
           end
           times = obj.Spikes.Times;
           sel_times = times(times > ops.cutout(1) & times < ops.cutout(2));
           
           figure('Color','w','Position',[100 500 1500 500]);
           if type == "hist"
               % histogram(sel_times,'BinWidth',ops.binning,'DisplayStyle','stairs')
               histogram(sel_times,ops.cutout(1):ops.binning:ops.cutout(2),'DisplayStyle','stairs')
           elseif type == "scatter"
                units = obj.Spikes.Units;
                sel_units = units(times > ops.cutout(1) & times < ops.cutout(2));
                scatter(sel_times,sel_units,0.5,'k.')

           else
               error('Unknown input')
           end
           % axis tight
           if ops.input == "raw"
               t_start = obj.Bursts.raw.T_start;
               t_end = obj.Bursts.raw.T_end;
           elseif ops.input == "processed"
               t_start = obj.Bursts.T_start;
               t_end = obj.Bursts.T_end;
           else
               warning("Unknown input type, using processed burst data")
           end
           if ~isempty(obj.Bursts.T_start)
               sel_t_start = t_start(t_start>ops.cutout(1) & t_start < ops.cutout(2));
               sel_t_end = t_end(t_end>ops.cutout(1) & t_end < ops.cutout(2));
               hold on
               xline(sel_t_start,'g--','LineWidth',1,'Label','Start','LabelHorizontalAlignment','center')
               xline(sel_t_end,'r--','LineWidth',1,'Label','End','LabelHorizontalAlignment','center')
           elseif isempty(obj.Bursts)
               warning('No burst data found. Run obj.detectBursts() first')
           end
           box off
           xlim(ops.cutout)
           % xlim((ops.cutout) /ops.binning)
           
       end
       
       function PlotCCG(obj,ccg_idx1,ccg_idx2, color)
           arguments
               obj MEArecording
               ccg_idx1
               ccg_idx2
               color = [0 0 0]
           end
           %            f = figure('Color','w');
           a1 = gca;%axes(f);
           dur = obj.Parameters.CCG.Duration * 2;
           [ccg,t] = CCG({obj.Units(ccg_idx1).SpikeTimes,obj.Units(ccg_idx2).SpikeTimes},[],'binSize',obj.Parameters.CCG.BinSize,...
               'duration',dur,'Fs',1/obj.RecordingInfo.SamplingRate);
           bar(a1,t,ccg(:,1,2),'k','BarWidth',1)
           if ccg_idx1 ~= ccg_idx2
               bar(a1,t,ccg(:,1,2),'FaceColor',color,'EdgeColor',color,'BarWidth',1)
           else
               bar(a1,t,ccg(:,1,1),'FaceColor',color,'EdgeColor',color,'BarWidth',1)
           end
           set(gca,'FontSize',7)
           
           xlabel('Interval (ms)')
%            xline(0,'r');
           a1.YAxis.Visible = 'off';
           box off
           axis tight
           xlims = xlim;
           xticks([xlims(1) 0 xlims(2)])
           x_max = dur/2*1000;
           xticklabels([-x_max 0 x_max])
       end
       
       function PlotCommunity(obj,method,matrix,flag)
           arguments
               obj MEArecording
               method string %"DDC","STTC","CCG"
               matrix string = "wu" %"wu" or "bu"
               flag logical = 1
           end
           %flag indicates if squares are drawnaround communities
           % M     = community_louvain(con_mat,0.0001,[],'negative_asym');
           if ~isfield(obj.Connectivity.(method), matrix)
               warning("No results available, running inference...")
               obj.inferConnectivity(method);
               disp("Inference done")
           end
           con_mat = obj.Connectivity.(method).(matrix);
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
       
       function data = PlotUnitClusterActivity(obj, time_cutout)
           ids = [obj.Units.ClusterID];
           clust_sz = histcounts(ids);
           N_clust = max(ids);
           clust_act = {};
           clust_ids = {};
           color_array = [];
           ytick_vals = [];
           for i = 1:N_clust
               id = length(clust_act) + 1;
               act = {obj.Units(ids == i).SpikeTimes};
               ytick_vals(i) = length(clust_act) + round(length(act)/2);
               clust_ids = [clust_ids arrayfun(@(x) ones(size(act{x}))* (x + length(clust_act)),1:length(act),'un',0)];
               clust_act = [clust_act act];
               color_array = [color_array ones(1,length(vertcat(act{:}))) * i];
           end
           spk_t = vertcat(clust_act{:});
           spk_id = vertcat(clust_ids{:});
           figure('Color','w'); 
           scatter(spk_t,spk_id,2,color_array,'filled')
           hold on
           yline(cumsum(clust_sz),'k-')
           xlim([time_cutout])
           yticks(ytick_vals)
           yticklabels(1:length(ytick_vals))
           ylim([0 sum(clust_sz)])
           colormap(othercolor('Mdarkrainbow',6))
           data.spk_t = spk_t;
           data.spk_id = spk_id;
           data.color_array = color_array;
           data.clust_sz = clust_sz;
           data.ytick_vals = ytick_vals;
       end
       
       function [colors, v] = PlotUnitTemplate(obj, unit_ids, scale_factor, wf_cutout, colors)
           arguments
              obj MEArecording
              unit_ids double
              scale_factor = 1.9
              wf_cutout = 1:60
              colors (:,3) double = othercolor('Set19',length(unit_ids));
           end
           
           template_matrix = obj.getTemplateMatrix();
           template_matrix = template_matrix(:,sum(template_matrix,[1,3],'omitnan')~=0,:); %Remove all zero paddings
           all_x = [];
           all_y = [];
%            figure('Color','w');
           hold on
           for u = 1:length(unit_ids)
               unit_template = squeeze(template_matrix(obj.Units(unit_ids(u)).TemplateID,wf_cutout,:));
               empty_channels = all(abs(unit_template)<0.02);
               clean_template = unit_template(:,~empty_channels);
               norm_template = clean_template./max(abs(clean_template),[],'all');
               coor = obj.getElectrodeCoordinates();
               x_ais{u} = coor(obj.Units(unit_ids(u)).ReferenceElectrode,1)/17.5;
               y_ais{u} = coor(obj.Units(unit_ids(u)).ReferenceElectrode,2)/17.5;
               
               clean_coor = coor(~empty_channels,:)/17.5;
               x{u} = clean_coor(:,1) + linspace(0,scale_factor,size(norm_template,1));
               y{u} = clean_coor(:,2) + (norm_template' * scale_factor);
               
               all_x = [x{u}; all_x];
               all_y = [y{u}; all_y];
           end
           x_els = min(mean(all_x,2)):2:max(mean(all_x,2));
           y_els = min(mean(all_y,2)):2:max(mean(all_y,2));
           [x_els_grid, y_els_grid] = ndgrid(x_els,y_els);
           scatter(x_els_grid, y_els_grid,70,'filled','Marker','square','MarkerFaceAlpha',0.1,'MarkerFaceColor','k',...
           'MarkerEdgeColor','k','MarkerEdgeAlpha',0.3)
           
           for u = 1:length(unit_ids)
               v{u} = viscircles([x_ais{u}+(scale_factor/2), y_ais{u}],1,'color',colors(u,:));
               plot(x{u}',y{u}','Color',colors(u,:),'LineWidth',1.5)
           end
       end
       
       function PlotUnitTraces(obj, raw_path, dataset, unit, cutout, spike_id, color, N_traces, filter_traces)
           arguments
               obj MEArecording
               raw_path string
               dataset string
               unit Unit
               cutout = 1000
               spike_id = 10
               color = 'r'
               N_traces = 4
               filter_traces = true
           end
           full_cutout = cutout * 2;
           spk_t = unit.SpikeTimes * unit.MEArecording.RecordingInfo.SamplingRate;
           tmp = unit.MEArecording.getTemplateMatrix();
           unit_tmp = squeeze(tmp(unit.TemplateID,:,:));
           max_tmp = max(unit_tmp);
           [~, sort_idx] = sort(max_tmp,'descend');
           
           raw_traces = h5read(raw_path,dataset,[spk_t(spike_id)-full_cutout, 1],[full_cutout*2, 1020]);
           unit_traces = raw_traces(:,sort_idx(1:N_traces));
           if filter_traces
               
               filt = designfilt('bandpassfir','FilterOrder',55,'CutoffFrequency1',300, 'CutoffFrequency2',4500,'SampleRate',10000);
               unit_traces = filtfilt(filt,double(unit_traces));
           end
           raw_offset = unit_traces + linspace(0,100,size(unit_traces,2));
           raw_offset = raw_offset(cutout:end-cutout,:);
           
           color_cutout = round(length(raw_offset)/2)-30:round(length(raw_offset)/2)+30;
           plot(1:length(raw_offset), raw_offset, 'k')
           hold on
           plot(color_cutout, raw_offset(color_cutout,:),'Color',color)
           axis tight
           axis off
           
       end
       
       function PlotActivityScan(obj, activity_scan)
           arguments
               obj MEArecording
               activity_scan string
           end
           well_info = h5info(activity_scan,"/wells/well001");
           data_store_names = cellfun(@string, {well_info.Links.Value});
           data_store_settings = h5read(activity_scan, [char(data_store_names(1)), '/settings/mapping']);
           data_store_data = h5read(activity_scan, [char(data_store_names(1)), '/spikes']);
           channel_spikes = groupcounts(data_store_data.channel);
       end
    end
end