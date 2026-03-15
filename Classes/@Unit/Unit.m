classdef Unit < handle
    properties
        MEArecording
        ReferenceWaveform
        ReferenceElectrode
        TemplateID
        SpikeTimes
        WaveformFeatures
        ActivityFeatures
        RegularityFeatures
        GraphFeatures
        Catch22
        ACG
    end

    properties (SetObservable = true)
        ClusterID
    end

    properties (Dependent)
        FullACG
    end

    properties (Access = private)
        RecordingDuration    % copied from MEArecording.RecordingInfo.Duration
        SamplingRate         % copied from MEArecording.RecordingInfo.SamplingRate
        UnitParams           % copied subset of MEArecording.Parameters (Regularity, Catch22, CCG, Analyses)
        CachedUnitID         % cached unit index in MEArecording.Units
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
            arguments
                mearec = []
                waveform = []
                reference_electrode = []
                spike_times = []
                waveform_features = []
            end
            if isempty(mearec)
                return % allow no-arg construction for array pre-allocation
            end
            unit.MEArecording = mearec;
            unit.ReferenceWaveform = waveform;
            unit.ReferenceElectrode = reference_electrode;
            unit.SpikeTimes = spike_times;
            unit.WaveformFeatures = struct2table(waveform_features);
            % Cache frequently accessed recording info to reduce parent traversal
            unit.RecordingDuration = mearec.RecordingInfo.Duration;
            unit.SamplingRate = mearec.RecordingInfo.SamplingRate;
            unit.UnitParams.Regularity = mearec.Parameters.Regularity;
            unit.UnitParams.Catch22 = mearec.Parameters.Catch22;
            unit.UnitParams.CCG = mearec.Parameters.CCG;
            unit.UnitParams.Analyses = mearec.Parameters.Analyses;
            unit.inferActivityFeatures();
        end

        function computeACG(unit, params)
        % COMPUTEACG  Compute and store the autocorrelogram.
        %   unit.computeACG()         — use default CCG params or extract from connectivity
        %   unit.computeACG(params)   — params struct with .BinSize, .Duration fields
            arguments
                unit Unit
                params struct = struct()
            end
            if isempty(fieldnames(params))
                try
                    acg = unit.MEArecording.Connectivity.CCG.CCGs(:, unit.unitID, unit.unitID);
                catch
                    [acg, ~] = CCG(unit.SpikeTimes, ones(size(unit.SpikeTimes)), ...
                        'binSize', unit.UnitParams.CCG.BinSize, ...
                        'duration', unit.UnitParams.CCG.Duration, ...
                        'Fs', 1/unit.SamplingRate);
                end
            else
                [acg, ~] = CCG(unit.SpikeTimes, ones(size(unit.SpikeTimes)), ...
                    'binSize', params.BinSize, ...
                    'duration', params.Duration, ...
                    'Fs', 1/unit.SamplingRate);
            end
            unit.ACG = double(acg);
        end

        function uid = get.unitID(unit)
            % Return cached ID; recompute only if cache is empty
            if isempty(unit.CachedUnitID) && ~isempty(unit.MEArecording)
                unit.CachedUnitID = find(unit.MEArecording.Units == unit, 1);
            end
            uid = unit.CachedUnitID;
        end

        function acg = get.FullACG(unit)
            acg = unit.MEArecording.Connectivity.FullCCG.CCGs(:, unit.unitID, unit.unitID);
        end
    end
end
