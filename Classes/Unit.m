classdef Unit < handle
    properties
        MEArecording
        ReferenceWaveform
        ReferenceElectrode
        TemplateID %Added
        SpikeTimes
        WaveformFeatures
        ActivityFeatures
        RegularityFeatures
        GraphFeatures
        Catch22
        ACG
    end

    properties(SetObservable = true)
        ClusterID
    end

    properties(Dependent)
        unitID
        FullACG
    end
    
    methods (Static)
                
        function feature_names = returnFeatureNames(feature_group)
            switch feature_group
                case "act"
                    feature_names = ["FiringRate","MeanInterSpikeInterval","VarianceInterSpikeInterval","CVInterSpikeInterval","PartialAutocorrelation"];
                case "reg"
                    feature_names = ["RegularityFrequency","RegularityMagnitude","RegularityFit"];
                case "c22"
                    feature_names = string(GetAllFeatureNames());
            end
        end
        
    end

    
    methods
        
        function unit = Unit(mearec, waveform, reference_electrode, spike_times, waveform_features)
            if nargin > 0
                unit.MEArecording = mearec;
                unit.ReferenceWaveform = waveform;
                unit.ReferenceElectrode = reference_electrode;
                unit.SpikeTimes = spike_times;
                unit.WaveformFeatures = struct2table(waveform_features);
                unit.inferActivityFeatures();
            end
        end

        function unit = set.ACG(unit,params)
            arguments
                unit Unit
                params struct = struct()
            end
            if isempty(params)
                try
                    acg = unit.MEArecording.Connectivity.CCG.CCGs(:,unit.unitID,unit.unitID);
                catch
                    [acg,t] = CCG(unit.SpikeTimes,ones(size(unit.SpikeTimes)),'binSize',unit.MEArecording.Parameters.CCG.BinSize,'duration',unit.MEArecording.Parameters.CCG.Duration,'Fs',1/unit.MEArecording.RecordingInfo.SamplingRate);
                end
            else
                [acg,t] = CCG(unit.SpikeTimes,ones(size(unit.SpikeTimes)),'binSize',params.BinSize,'duration',params.Duration,'Fs',1/unit.MEArecording.RecordingInfo.SamplingRate);
            end
            unit.ACG = double(acg);
        end

        function unit_id = get.unitID(unit)
            unit_id = find(unit.MEArecording.Units == unit);
        end

        % function acg = get.ACG(unit)
        % 
        %     % acg = unit.MEArecording.Connectivity.CCG.CCGs(:,unit.unitID,unit.unitID);
        % end

        function acg = get.FullACG(unit)
            acg = unit.MEArecording.Connectivity.FullCCG.CCGs(:,unit.unitID,unit.unitID);
        end
        
        function [eCCG,iCCG] = getCCGconnections(unit)
            % unit.MEArecording.Connectivity.C
            % e_idx = 
        end
        
        function inferActivityFeatures(unit)
            if length(unit.SpikeTimes) > 2
                isi = diff(unit.SpikeTimes);
                act_feat.FiringRate = length(unit.SpikeTimes)/unit.MEArecording.RecordingInfo.Duration;
                act_feat.MeanInterSpikeInterval = mean(isi);
                act_feat.VarianceInterSpikeInterval = var(isi);
                act_feat.CVInterSpikeInterval = std(isi)/mean(isi);
                pacf = parcorr(isi,1);
                act_feat.PartialAutocorrelation = pacf(2);
                unit.ActivityFeatures = struct2table(act_feat);
                if unit.MEArecording.Parameters.Analyses.Regularity
                    unit.RegularityFeatures = unit.getRegularity();
                end
                if unit.MEArecording.Parameters.Analyses.Catch22
                    catch_22_table = unit.MEArecording.run_catch_22(unit.SpikeTimes);
                    catch_22_table.Properties.VariableNames = "SC_" + string(catch_22_table.Properties.VariableNames);
                    unit.Catch22 = catch_22_table;
                    
                end
            end
        end
        
        function regularity_table = getRegularity(unit)
           binned_activity = histcounts(unit.SpikeTimes,'BinLimits',[0 unit.MEArecording.RecordingInfo.Duration],'BinWidth',unit.MEArecording.Parameters.Regularity.Binning);
           norm_activity = binned_activity/max(binned_activity);
           NFFT = length(norm_activity);
           F = (0 : 1/NFFT : 1/2-1/NFFT)*(1/unit.MEArecording.Parameters.Regularity.Binning);
           TEMP = fft(norm_activity,NFFT);
           TEMP(1) = 0;
           freq_domain = abs(TEMP(1:NFFT/2));
           [mag,freq_idx] = max(freq_domain);
           reg_feat.RegularityFrequency = F(freq_idx);
           reg_feat.RegularityMagnitude = mag;
           
           norm_freq_domain = freq_domain/max(freq_domain);
           l = [];
           p = [];
           for i = 1:unit.MEArecording.Parameters.Regularity.N_peaks
               [pks,locs,~,~] = findpeaks(norm_freq_domain,'NPeaks',1,'SortStr','descend');
               norm_freq_domain = norm_freq_domain(locs:end);
               l = [l; locs]; %#ok<AGROW>
               p = [p; pks]; %#ok<AGROW>
               if isempty(pks) || length(norm_freq_domain) < 10
                   break
               end
           end
           log_p = log10(p)-min(log10(p)); % Make data positive to be able to fit
           try
               f = fit(cumsum(l),log_p,'exp1');
               reg_feat.RegularityFit = f.b;
           catch
               reg_feat.RegularityFit = NaN;
           end
           regularity_table = struct2table(reg_feat);
       end
       
    end
end