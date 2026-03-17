classdef Unit < handle

    properties (Constant, Hidden)
        ClassVersion = 1  % Increment on each schema change
    end

    properties
        MEArecording
        ReferenceWaveform
        ReferenceElectrode
        TemplateID
        StableID string     % Deterministic ID: hash of (ChipID, PlatingDate, TemplateID, ReferenceElectrode)
        SpikeTimes
        WaveformFeatures
        ActivityFeatures
        RegularityFeatures
        GraphFeatures
        Catch22
        ACG
        BombcellMetrics   % Table: per-unit bombcell quality metrics (BC_ prefixed columns)
        BombcellType      % Scalar: bombcell classification (0=noise, 1=good, 2=MUA)
        KSLabel string    % Kilosort label from cluster_KSLabel.tsv ("good" or "mua")
    end

    properties (SetObservable = true)
        ClusterID
    end

    properties (Dependent)
        FullACG
        unitID           % Index of this unit in MEArecording.Units (recomputed each access)
    end

    properties (Access = private)
        RecordingDuration    % copied from MEArecording.RecordingInfo.Duration
        SamplingRate         % copied from MEArecording.RecordingInfo.SamplingRate
        UnitParams           % copied subset of MEArecording.Parameters (Regularity, Catch22, CCG, Analyses)
        CachedUnitID         % cached unit index in MEArecording.Units
    end

    properties (Access = private, Transient)
        FeatureCache         % containers.Map keyed by "featureName_paramHash" — not serialized
    end

    methods (Static)

        function feature_names = returnFeatureNames(feature_group)
            switch feature_group
                case "act"
                    feature_names = ["FiringRate","MeanInterSpikeInterval","VarianceInterSpikeInterval","CVInterSpikeInterval", ...
                        "CV2InterSpikeInterval","LocalVariation","RevisedLocalVariation","FanoFactor","PartialAutocorrelation"];
                case "reg"
                    feature_names = ["RegularityFrequency","RegularityMagnitude","RegularityFit"];
                case "c22"
                    feature_names = string(GetAllFeatureNames());
                case "bc"
                    feature_names = ["BC_nPeaks","BC_nTroughs","BC_waveformDuration", ...
                        "BC_mainPeakToTroughRatio","BC_peak1ToPeak2Ratio", ...
                        "BC_scndPeakToTroughRatio","BC_troughToPeak2Ratio", ...
                        "BC_waveformBaselineFlatness","BC_spatialDecaySlope", ...
                        "BC_fractionRPVs","BC_presenceRatio","BC_percSpikesMissing"];
            end
        end

        function obj = loadobj(s)
        % LOADOBJ  Reconstruct Unit from saved data.
            if isstruct(s)
                obj = Unit();
                fnames = fieldnames(s);
                for i = 1:length(fnames)
                    try
                        obj.(fnames{i}) = s.(fnames{i});
                    catch
                        % Skip fields that no longer exist or are read-only
                    end
                end
            else
                obj = s;
            end
            % Backfill StableID for legacy objects
            if isempty(obj.StableID) && ~isempty(obj.MEArecording)
                obj.generateStableID();
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
            unit.FeatureCache = containers.Map();
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
            unit.UnitParams.QC.RefractoryPeriod = mearec.Parameters.QC.RefractoryPeriod;
            unit.generateStableID();
            unit.inferActivityFeatures();
        end

        function ensureCache(unit)
        % ENSURECACHE  Lazily initialize FeatureCache (empty after load since Transient).
            if isempty(unit.FeatureCache)
                unit.FeatureCache = containers.Map();
            end
        end

        function generateStableID(unit)
        % GENERATESTABLEID  Create deterministic identifier from recording metadata + unit identity.
        %   ID is a hex string derived from (ChipID, PlatingDate, TemplateID, ReferenceElectrode).
        %   Deterministic: same inputs always produce the same ID.
            if isempty(unit.MEArecording)
                return
            end
            md = unit.MEArecording.Metadata;
            chip = "unknown";
            plating = "unknown";
            if isfield(md, 'ChipID'),      chip    = string(md.ChipID);      end
            if isfield(md, 'PlatingDate'),  plating = string(md.PlatingDate); end
            key = sprintf('%s_%s_%d_%d', chip, plating, unit.TemplateID, unit.ReferenceElectrode);
            unit.StableID = string(dec2hex(abs(java.lang.String(key).hashCode()), 8));
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
                has_ccg = ~isempty(unit.MEArecording) ...
                    && isfield(unit.MEArecording.Connectivity, 'CCG') ...
                    && isfield(unit.MEArecording.Connectivity.CCG, 'CCGs') ...
                    && ~isempty(unit.unitID);
                if has_ccg
                    acg = unit.MEArecording.Connectivity.CCG.CCGs(:, unit.unitID, unit.unitID);
                else
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
            % Always recompute to avoid stale indices after unit filtering/removal
            if ~isempty(unit.MEArecording) && ~isempty(unit.MEArecording.Units)
                uid = find(unit.MEArecording.Units == unit, 1);
            else
                uid = [];
            end
        end

        function acg = get.FullACG(unit)
            if isempty(unit.MEArecording) || ~isfield(unit.MEArecording.Connectivity, 'FullCCG')
                acg = [];
                return
            end
            acg = unit.MEArecording.Connectivity.FullCCG.CCGs(:, unit.unitID, unit.unitID);
        end

        function s = saveobj(unit)
        % SAVEOBJ  Serialize Unit for save(), breaking the circular MEArecording reference.
            s = struct();
            mc = metaclass(unit);
            props = mc.PropertyList;
            for i = 1:length(props)
                p = props(i);
                if p.Dependent || p.Transient || ~strcmp(p.GetAccess, 'public')
                    continue
                end
                if p.Name == "MEArecording"
                    continue  % Break circular reference — restored by MEArecording.loadobj
                end
                s.(p.Name) = unit.(p.Name);
            end
        end
    end
end
