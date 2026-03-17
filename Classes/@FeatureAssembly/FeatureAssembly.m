classdef FeatureAssembly
% FEATUREASSEMBLY  Static methods for assembling feature tables from objects.
%
%   Consolidates feature assembly logic previously scattered across
%   MEArecording.getUnitFeatures, MEArecording.getRecordingFeatures,
%   RecordingGroup.aggregateCultureFeatureTables, and RecordingGroup.alignWaveforms.
%
%   All methods are Static — no instance needed.
%
%   USAGE:
%     feat_table = FeatureAssembly.unitFeatures(unit_array, ["ActivityFeatures","WaveformFeatures"])
%     feat_table = FeatureAssembly.recordingFeatures(recordings, "all", "all", false)
%     [feat_table, cultures] = FeatureAssembly.cultureFeatures(rg, "Recording", "DIV", [7 14 21 28])

    methods (Static)

        function feat_table = unitFeatures(obj, unit_features)
        % UNITFEATURES  Build a table of per-unit features from recordings or unit array.
        %
        % INPUTS:
        %   obj           - MEArecording array or Unit array
        %   unit_features - string array of feature groups to include
        %                   (e.g. "ActivityFeatures", "ReferenceWaveform", "ACG", "all")
            arguments
                obj
                unit_features string = ["ReferenceWaveform", "ActivityFeatures"]
            end
            if isscalar(unit_features) && unit_features == "all"
                unit_features = ["ActivityFeatures", "WaveformFeatures", "RegularityFeatures", "Catch22", "GraphFeatures"];
            end
            if isa(obj, 'MEArecording')
                unit_array = [obj.Units];
                % Auto-include BombcellMetrics if available
                if ~isempty(unit_array) && ~isempty(unit_array(1).BombcellMetrics) && ~any(unit_features == "BombcellMetrics")
                    unit_features = [unit_features, "BombcellMetrics"];
                end
            elseif isa(obj, 'Unit')
                unit_array = obj;
            else
                error('FeatureAssembly:unitFeatures', 'Input must be MEArecording array or Unit array.');
            end
            feature_tables = {};
            for f = 1:length(unit_features)
                if unit_features(f) == "ReferenceWaveform"
                    aligned_wf = FeatureAssembly.alignWaveforms(unit_array);
                    var_names = "Waveform" + (1:size(aligned_wf, 1));
                    feature_tables{f} = array2table(aligned_wf', "VariableNames", var_names);
                elseif ismember(unit_features(f), ["ACG", "FullACG"])
                    acgs = [unit_array.(unit_features(f))];
                    norm_acgs = acgs ./ max(acgs);
                    norm_acgs(isnan(norm_acgs)) = 0;
                    var_names = unit_features(f) + (1:size(norm_acgs, 1));
                    feature_tables{f} = array2table(norm_acgs', "VariableNames", var_names);
                else
                    feature_tables{f} = vertcat(unit_array.(unit_features(f)));
                end
            end
            feat_table = [feature_tables{:}];
        end

        function feat_table = recordingFeatures(obj, network_features, unit_features, useClustered)
        % RECORDINGFEATURES  Build a table of per-recording features.
        %
        % INPUTS:
        %   obj              - MEArecording array
        %   network_features - "all" or string array of network feature names, or []
        %   unit_features    - "all" or string array of unit feature names, or []
        %   useClustered     - logical, use clustered features
            arguments
                obj MEArecording
                network_features string = "all"
                unit_features string = "all"
                useClustered logical = false
            end

            if ~isempty(network_features)
                fname_array = arrayfun(@(x) fieldnames(x.NetworkFeatures), obj, 'un', 0);
                [gc, fnames] = groupcounts(vertcat(fname_array{:}));
                rm_field = string(fnames{gc ~= length(fname_array)});
                fnames = fnames(gc == length(fname_array));
                for i = 1:length(obj)
                    network_struct = obj(i).NetworkFeatures;
                    if isfield(network_struct, rm_field)
                        network_struct = rmfield(network_struct, rm_field);
                    end
                    network_array(i) = network_struct;
                end
                fnames = fieldnames(network_array);
                if isscalar(network_features) && network_features == "all"
                    nw_idx = 1:length(fnames);
                else
                    nw_idx = find(contains(fnames, network_features));
                end
                network_cell = squeeze(struct2cell(network_array));
                network_cell = reshape(network_cell, [], length(obj));
                for i = 1:size(network_cell, 2)
                    network_table(i,:) = [network_cell{nw_idx, i}];
                end
            else
                network_table = table();
            end

            if ~isempty(unit_features)
                if useClustered
                    unit_table = concatenateClusteredFeatures(obj, 0, unit_features);
                else
                    fnames = fieldnames([obj.UnitFeatures]);
                    if isscalar(unit_features) && unit_features == "all"
                        feat_idx = 1:length(fnames);
                    else
                        feat_idx = find(contains(fnames, unit_features));
                        if length(feat_idx) ~= length(unit_features)
                            error(join('Could not find: ' + unit_features(~contains(unit_features, fnames))))
                        end
                    end
                    unit_cell = squeeze(struct2cell([obj.UnitFeatures]));
                    unit_cell = reshape(unit_cell, [], length(obj));
                    for i = 1:size(unit_cell, 2)
                        unit_table(i,:) = [unit_cell{feat_idx, i}];
                    end
                end
            else
                unit_table = table();
            end

            feat_table = [network_table unit_table];
        end

        function [feature_table, culture_array, values] = cultureFeatures(rg, level, metadata_field, values, tolerance, network_features, unit_features, normalization, useClustered)
        % CULTUREFEATURES  Aggregate features across culture timepoints.
        %
        % Wrapper that delegates to the existing aggregateCultureFeatureTables
        % logic. Will be refactored further in future iterations.
            arguments
                rg RecordingGroup
                level string = "Recording"
                metadata_field string = "DIV"
                values {isnumeric} = [7, 14, 21, 28]
                tolerance {isnumeric} = 1
                network_features string = "all"
                unit_features string = "all"
                normalization string = "baseline"
                useClustered logical = false
            end
            [feature_table, culture_array, values] = aggregateCultureFeatureTables(rg, level, metadata_field, values, tolerance, network_features, unit_features, normalization, useClustered);
        end

        function aligned_wf = alignWaveforms(unit_array, target_sr)
        % ALIGNWAVEFORMS  Interpolate and trough-align reference waveforms.
        %
        % INPUTS:
        %   unit_array - (1 x N) Unit array
        %   target_sr  - target sampling rate in Hz after interpolation
            arguments
                unit_array Unit
                target_sr  (1,1) double = 0
            end

            actual_sr = unit_array(1).MEArecording.RecordingInfo.SamplingRate;
            if target_sr == 0
                target_sr = actual_sr * 10;
            end
            interp_factor_raw = target_sr / actual_sr;
            interp_factor     = round(interp_factor_raw);
            if abs(interp_factor_raw - interp_factor) > 1e-9
                warning('FeatureAssembly:alignWaveforms', ...
                    'target_sr (%g Hz) is not an integer multiple of the recording rate (%g Hz). ' ...
                    'Rounding interpolation factor from %.6f to %i (effective rate: %g Hz).', ...
                    target_sr, actual_sr, interp_factor_raw, interp_factor, actual_sr * interp_factor);
            end

            ref_wf = [unit_array.ReferenceWaveform];
            ref_wf = ref_wf(sum(ref_wf, 2) ~= 0, :);
            [~, i] = min(ref_wf, [], 1);
            peak_idx = mean(i);
            max_offset = round(peak_idx / 2);
            x = max_offset:size(ref_wf, 1) + max_offset - 1;
            buffer_wf = size(ref_wf, 1) + 2 * max_offset;
            xq = linspace(1, buffer_wf, buffer_wf * interp_factor);

            interp_wf = interp1(x, ref_wf, xq, "makima");
            interp_wf = interp_wf((max_offset * interp_factor) + 1:((buffer_wf - max_offset - 1) * interp_factor), :);
            interp_wf = interp_wf ./ max(abs(interp_wf));
            rm_idx = find(ceil(abs(peak_idx - i)) >= max_offset);
            [~, minIndices] = min(interp_wf);
            shifts = round(peak_idx * interp_factor) - minIndices;

            aligned_wf = zeros(size(interp_wf));
            for j = 1:size(interp_wf, 2)
                if ismember(j, rm_idx)
                    aligned_wf(:, j) = interp_wf(:, j);
                else
                    aligned_wf(:, j) = circshift(interp_wf(:, j), [shifts(j), 0]);
                end
            end
            shifts(rm_idx) = [];
            aligned_wf = aligned_wf((5 * interp_factor):(end - 5 * interp_factor), :);
        end

    end
end
