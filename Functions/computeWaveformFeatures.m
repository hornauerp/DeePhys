function waveform_features = computeWaveformFeatures(norm_wf_matrix, sampling_rate)
% COMPUTEWAVEFORMFEATURES  Compute waveform shape features from normalised waveforms.
%
%   waveform_features = computeWaveformFeatures(norm_wf_matrix, sampling_rate)
%
% INPUTS:
%   norm_wf_matrix - (N_samples x N_units) double, normalised peak-channel waveforms
%   sampling_rate  - (1 x 1) double, sampling rate in Hz
%
% OUTPUTS:
%   waveform_features - (1 x N_units) cell array of structs, each with fields:
%       AUC_peak_1, AUC_trough, AUC_peak_2, Rise, Decay, HalfWidth,
%       Asymmetry, T2Pdelay, T2Pratio

    arguments
        norm_wf_matrix (:,:) double
        sampling_rate  (1,1) double
    end

    interpolation_factor = 10;
    ms_conversion = 10 * interpolation_factor;
    x = 1:size(norm_wf_matrix,1);
    xq = 1:(1/interpolation_factor):size(norm_wf_matrix,1);
    interp_wf_matrix = double(interp1(x, norm_wf_matrix, xq, 'pchip'));
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    n_units = size(interp_wf_matrix, 2);

    unit_zero_crossings = cell(1, n_units);
    for u = 1:n_units
        zc_idx = zci(interp_wf_matrix(:,u));
        unit_zero_crossings{u} = zc_idx(:)';
    end

    [unit_trough_value, unit_trough_idx] = min(interp_wf_matrix);
    waveform_features = cell(1, n_units);

    for u = 1:n_units
        unit_features = struct();

        p1_start = max(1, unit_trough_idx(u) - ms_conversion);
        p1_end   = unit_trough_idx(u);
        p2_start = unit_trough_idx(u);
        p2_end   = min(size(interp_wf_matrix,1), unit_trough_idx(u) + ms_conversion);

        peak_1_cutout = interp_wf_matrix(p1_start:p1_end, u);
        peak_2_cutout = interp_wf_matrix(p2_start:p2_end, u);

        [peak_1_value, ~] = max(peak_1_cutout);
        [peak_2_value, peak_2_rel_idx] = max(peak_2_cutout);
        peak_2_abs_idx = peak_2_rel_idx + unit_trough_idx(u) - 1;

        denom_asym  = peak_2_value + peak_1_value;
        asymmetry_u = (peak_2_value - peak_1_value) / denom_asym;
        if denom_asym == 0, asymmetry_u = NaN; end
        t2pdelay_u  = (peak_2_abs_idx - unit_trough_idx(u)) / (sampling_rate * interpolation_factor / 1000);
        if peak_2_value == 0
            t2pratio_u = NaN;
        else
            t2pratio_u = abs(unit_trough_value(u) / peak_2_value);
        end

        zc = [1 unit_zero_crossings{u} length(xq)];
        [~,zc_pre_trough] = max(zc(zc<unit_trough_idx(u)));
        zc_pre_trough = max([zc_pre_trough 2]);
        unit_features.AUC_peak_1 = trapz(interp_wf_matrix(zc(zc_pre_trough - 1):zc(zc_pre_trough), u));
        if zc_pre_trough + 1 > length(zc)
            unit_features.AUC_trough = NaN;
            unit_features.AUC_peak_2 = NaN;
        else
            unit_features.AUC_trough = abs(trapz(interp_wf_matrix(zc(zc_pre_trough):zc(zc_pre_trough + 1), u)));
            if length(zc(zc_pre_trough:end)) > 3
                unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):zc(zc_pre_trough + 2), u));
            else
                unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):size(interp_wf_matrix,1), u));
            end
        end

        rise_cutout = interp_wf_matrix(unit_trough_idx(u):peak_2_abs_idx, u);
        padded_rise = [ones(100,1)*rise_cutout(1); rise_cutout; ones(100,1)*rise_cutout(end)];
        unit_features.Rise = mean(slewrate(padded_rise, interpolation_factor));
        decay_cutout = interp_wf_matrix(peak_2_abs_idx:end, u);
        decay_cutout = decay_cutout(1:find(decay_cutout==min(decay_cutout), 1, 'first'));
        padded_decay = [ones(100,1)*decay_cutout(1); decay_cutout; ones(100,1)*decay_cutout(end)];
        unit_features.Decay = mean(slewrate(padded_decay, interpolation_factor));

        half_amp = unit_trough_value(u) / 2;
        wf = interp_wf_matrix(:, u);
        above_half = wf < half_amp;
        if any(above_half)
            half_start = find(above_half, 1, 'first');
            half_end = find(above_half, 1, 'last');
            unit_features.HalfWidth = (half_end - half_start) / (sampling_rate * interpolation_factor / 1000);
        else
            unit_features.HalfWidth = NaN;
        end

        unit_features.Asymmetry = asymmetry_u;
        unit_features.T2Pdelay = t2pdelay_u;
        unit_features.T2Pratio = t2pratio_u;
        feature_names = fieldnames(unit_features);
        for f = 1:length(feature_names)
           if isempty(unit_features.(feature_names{f}))
               unit_features.(feature_names{f}) = NaN;
           end
        end
        waveform_features{u} = unit_features;
    end
end
