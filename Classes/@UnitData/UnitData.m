classdef UnitData
% UNITDATA  Per-unit raw data container. Value class, no back-reference to parent.
%
% UnitData holds the identity, raw waveform, and spike times for one sorted unit.
% It references its parent recording by RecordingID (a string key), not an object
% handle — eliminating the circular reference that required saveobj/loadobj hacks
% in the old Unit ↔ MEArecording design.
%
% USAGE:
%   ud = UnitData(template_id, ref_electrode, spike_times, ref_waveform, ...
%                 recording_id, sampling_rate, duration);
%   ud = UnitData.fromLegacyUnit(unit_obj, recording_id);
%   ud = UnitData.fromStruct(s);

    properties
        UnitID              string  % Stable deterministic ID (hash of ChipID+PlatingDate+TemplateID+RefElectrode)
        RecordingID         string  % Links to parent SpikeData.RecordingID (string, not a handle)
        TemplateID          double  % 1-indexed template number from Kilosort
        ReferenceElectrode  double  % Channel index of peak-amplitude electrode
        SpikeTimes          double  % (N_spikes x 1) spike times in seconds
        ReferenceWaveform   double  % (N_samples x 1) peak-channel waveform (raw, not interpolated)
        SamplingRate        double  % Copied from parent recording
        RecordingDuration   double  % Copied from parent recording (seconds)
        KSLabel             string  % Kilosort label ("good" / "mua"), "" if not available
        BombcellType        double  % BC classification: 0=noise, 1=good, 2=MUA, NaN=not run
    end

    methods

        function ud = UnitData(template_id, ref_electrode, spike_times, ref_waveform, ...
                               recording_id, sampling_rate, duration)
        % UNITDATA  Construct from raw inputs. All args optional for array pre-alloc.
            arguments
                template_id     double  = []
                ref_electrode   double  = []
                spike_times     double  = []
                ref_waveform    double  = []
                recording_id    string  = ""
                sampling_rate   double  = []
                duration        double  = []
            end
            if isempty(template_id)
                return   % empty construction for array pre-allocation
            end
            ud.TemplateID         = template_id;
            ud.ReferenceElectrode = ref_electrode;
            ud.SpikeTimes         = spike_times;
            ud.ReferenceWaveform  = ref_waveform;
            ud.RecordingID        = recording_id;
            ud.SamplingRate       = sampling_rate;
            ud.RecordingDuration  = duration;
            ud.KSLabel            = "";
            ud.BombcellType       = NaN;
        end

    end

    methods (Static)

        function ud = fromStruct(s)
        % FROMSTRUCT  Reconstruct UnitData from a plain struct (e.g. loaded from .mat).
            ud = UnitData();
            fns = fieldnames(s);
            for i = 1:numel(fns)
                try
                    ud.(fns{i}) = s.(fns{i});
                catch
                end
            end
        end

        function ud = fromLegacyUnit(unit_obj, recording_id)
        % FROMLEGACYUNIT  Extract UnitData from an old Unit handle object.
        %   Used by RecordingProcessor.fromLegacyMat.
            arguments
                unit_obj     Unit
                recording_id string
            end
            ud = UnitData();
            ud.UnitID             = unit_obj.StableID;
            ud.RecordingID        = recording_id;
            ud.TemplateID         = unit_obj.TemplateID;
            ud.ReferenceElectrode = unit_obj.ReferenceElectrode;
            ud.SpikeTimes         = unit_obj.SpikeTimes;
            ud.ReferenceWaveform  = unit_obj.ReferenceWaveform;

            % Use private-cached copies to avoid re-traversing MEArecording handle
            ud.SamplingRate       = unit_obj.SamplingRate;
            ud.RecordingDuration  = unit_obj.RecordingDuration;
            ud.KSLabel            = unit_obj.KSLabel;

            if ~isempty(unit_obj.BombcellType)
                ud.BombcellType = unit_obj.BombcellType;
            else
                ud.BombcellType = NaN;
            end
        end

        function ud_array = fromLegacyUnitArray(unit_array, recording_id)
        % FROMLEGACYUNITARRAY  Convert a Unit array to a UnitData array.
            ud_array = arrayfun(@(u) UnitData.fromLegacyUnit(u, recording_id), unit_array);
        end

    end
end
