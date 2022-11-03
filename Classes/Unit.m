classdef Unit < handle
    properties
        MEArecording
        ReferenceWaveform
        ReferenceElectrode
        SpikeTimes
        WaveformFeatures
        ActivityFeatures
    end
    
    methods
        function unit = Unit(mearec, waveform, reference_electrode, spike_times, waveform_features)
            if nargin > 0
                unit.MEArecording = mearec;
                unit.ReferenceWaveform = waveform;
                unit.ReferenceElectrode = reference_electrode;
                unit.SpikeTimes = spike_times;
                unit.WaveformFeatures = waveform_features;
                unit.inferActivityFeatures();
            end
        end
        
        function inferActivityFeatures(unit)
            for u = unit
                isi = diff(u.SpikeTimes);
                u.ActivityFeatures.MeanInterSpikeInterval = mean(isi);
                u.ActivityFeatures.VarianceInterSpikeInterval = std(isi);
                u.ActivityFeatures.CVInterSpikeInterval = std(isi)/mean(isi);
                pacf = parcorr(isi,1);
                u.ActivityFeatures.PartialAutocorrelation = pacf(2);
                u.getRegularity();
            end
        end
        
        function getRegularity(unit)
           binned_activity = histcounts(unit.SpikeTimes,'BinLimits',[0 unit.MEArecording.RecordingInfo.Duration],'BinWidth',unit.MEArecording.Parameters.Regularity.Binning);
           norm_activity = binned_activity/max(binned_activity);
           NFFT = length(norm_activity);
           F = (0 : 1/NFFT : 1/2-1/NFFT)*(1/unit.MEArecording.Parameters.Regularity.Binning);
           TEMP = fft(norm_activity,NFFT);
           TEMP(1) = 0;
           freq_domain = abs(TEMP(1:NFFT/2));
           [mag,freq_idx] = max(freq_domain);
           unit.ActivityFeatures.RegularityFrequency = F(freq_idx);
           unit.ActivityFeatures.RegularityMagnitude = mag;
           
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
               unit.ActivityFeatures.RegularityFit = f.b;
           catch
               unit.ActivityFeatures.RegularityFit = NaN;
           end
       end
       
    end
end