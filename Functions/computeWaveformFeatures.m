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
    ms_conversion = round(sampling_rate / 1000 * interpolation_factor);   % interpolated samples per ms (was hard-coded for 10 kHz)
    dt_ms = 1 / ms_conversion;   % AUC_* in amplitude x ms, the same at every sampling rate
    x = 1:size(norm_wf_matrix,1);
    xq = 1:(1/interpolation_factor):size(norm_wf_matrix,1);
    % Cubic spline: trough and peak can fall between samples. pchip is shape-preserving and
    % keeps every extremum on a sample, so the timing features stayed on the 1/fs grid.
    interp_wf_matrix = double(interp1(x(:), norm_wf_matrix, xq(:), 'spline'));
    if isvector(norm_wf_matrix), interp_wf_matrix = interp_wf_matrix(:); end   % interp1 returns a row for a single unit
    zci = @(v) find(v(1:end-1) .* v(2:end) <= 0);
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
        % Repolarisation peak searched up to 2 ms after the trough: broad units peak later
        % than 1 ms, and a 1 ms cap returned a point on the rising phase.
        p2_end   = min(size(interp_wf_matrix,1), unit_trough_idx(u) + 2 * ms_conversion);

        peak_1_cutout = interp_wf_matrix(p1_start:p1_end, u);
        peak_2_cutout = interp_wf_matrix(p2_start:p2_end, u);

        [peak_1_value, ~] = max(peak_1_cutout);
        [peak_2_value, peak_2_rel_idx] = max(peak_2_cutout);
        peak_2_abs_idx = peak_2_rel_idx + unit_trough_idx(u) - 1;

        denom_asym   = peak_2_value + peak_1_value;
        max_peak_abs = max(abs([peak_1_value, peak_2_value]));
        if max_peak_abs < 1e-9 || abs(denom_asym) < 1e-6 * max_peak_abs
            asymmetry_u = NaN;
        else
            asymmetry_u = (peak_2_value - peak_1_value) / denom_asym;
        end
        t2pdelay_u  = (peak_2_abs_idx - unit_trough_idx(u)) / (sampling_rate * interpolation_factor / 1000);
        if peak_2_value == 0
            t2pratio_u = NaN;
        else
            t2pratio_u = abs(unit_trough_value(u) / peak_2_value);
        end

        zc = [1 unit_zero_crossings{u} length(xq)];
        [~,zc_pre_trough] = max(zc(zc<unit_trough_idx(u)));
        zc_pre_trough = max([zc_pre_trough 2]);
        unit_features.AUC_peak_1 = dt_ms * trapz(interp_wf_matrix(zc(zc_pre_trough - 1):zc(zc_pre_trough), u));
        if zc_pre_trough + 1 > length(zc)
            unit_features.AUC_trough = NaN;
            unit_features.AUC_peak_2 = NaN;
        else
            unit_features.AUC_trough = dt_ms * abs(trapz(interp_wf_matrix(zc(zc_pre_trough):zc(zc_pre_trough + 1), u)));
            if length(zc(zc_pre_trough:end)) > 3
                unit_features.AUC_peak_2 = dt_ms * trapz(interp_wf_matrix(zc(zc_pre_trough + 1):zc(zc_pre_trough + 2), u));
            else
                unit_features.AUC_peak_2 = dt_ms * trapz(interp_wf_matrix(zc(zc_pre_trough + 1):size(interp_wf_matrix,1), u));
            end
        end

        rise_cutout = interp_wf_matrix(unit_trough_idx(u):peak_2_abs_idx, u);
        if numel(rise_cutout) > 1
            unit_features.Rise = mean(slewrate(rise_cutout, ms_conversion));   % amplitude per ms, independent of sampling rate
        else
            unit_features.Rise = NaN;
        end
        decay_cutout = interp_wf_matrix(peak_2_abs_idx:end, u);
        decay_cutout = decay_cutout(1:find(decay_cutout==min(decay_cutout), 1, 'first'));
        if numel(decay_cutout) > 1
            unit_features.Decay = mean(slewrate(decay_cutout, ms_conversion));
        else
            unit_features.Decay = NaN;
        end

        % Width of the trough at half its amplitude: the contiguous run below half-amplitude
        % that contains the trough (not the first and last such samples anywhere, which a
        % later negative deflection would join)
        half_amp = unit_trough_value(u) / 2;
        wf = interp_wf_matrix(:, u);
        above_half = wf < half_amp;
        ti = unit_trough_idx(u);
        if above_half(ti)
            half_start = find(~above_half(1:ti), 1, 'last') + 1;
            if isempty(half_start), half_start = 1; end
            half_end = ti - 1 + find(~above_half(ti:end), 1, 'first') - 1;
            if isempty(half_end), half_end = numel(wf); end
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
