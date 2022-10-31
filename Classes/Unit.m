classdef Unit < handle
    properties
        MEArecording
        ReferenceWaveform
        SpikeTimes
        WaveformFeatures
        ActivityFeatures
    end
    
    methods
        function unit = Unit(mearec, waveform, spike_times, waveform_features)
            if nargin > 0
                unit.MEArecording = mearec;
                unit.ReferenceWaveform = waveform;
                unit.SpikeTimes = spike_times;
                unit.WaveformFeatures = waveform_features;
                unit.inferActivityFeatures();
                
            end
        end
        
        function inferActivityFeatures(obj)
            
        end
    end
end