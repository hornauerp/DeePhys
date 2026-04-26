function aligned_wf = alignWaveformsFromData(ref_wf_matrix, sampling_rate, target_sr)
% ALIGNWAVEFORMSFROMDATA  Interpolate and trough-align reference waveforms.
%
%   aligned_wf = alignWaveformsFromData(ref_wf_matrix, sampling_rate)
%   aligned_wf = alignWaveformsFromData(ref_wf_matrix, sampling_rate, target_sr)
%
% INPUTS:
%   ref_wf_matrix - (N_samples x N_units) double, normalised peak-channel waveforms
%   sampling_rate - (1 x 1) double, recording sampling rate in Hz
%   target_sr     - (1 x 1) double, target sampling rate after interpolation
%                   (default: sampling_rate * 10)
%
% OUTPUTS:
%   aligned_wf - (M_samples x N_units) double, interpolated and trough-aligned waveforms

    arguments
        ref_wf_matrix (:,:) double
        sampling_rate (1,1) double
        target_sr     (1,1) double = 0
    end

    if target_sr == 0
        target_sr = sampling_rate * 10;
    end
    interp_factor_raw = target_sr / sampling_rate;
    interp_factor     = round(interp_factor_raw);
    if abs(interp_factor_raw - interp_factor) > 1e-9
        warning('alignWaveformsFromData:interpFactor', ...
            ['target_sr (%g Hz) is not an integer multiple of the recording rate (%g Hz). ' ...
            'Rounding interpolation factor from %.6f to %i (effective rate: %g Hz).'], ...
            target_sr, sampling_rate, interp_factor_raw, interp_factor, sampling_rate * interp_factor);
    end

    ref_wf = ref_wf_matrix(sum(ref_wf_matrix, 2) ~= 0, :);
    [~, i] = min(ref_wf, [], 1);
    peak_idx = mean(i);
    max_offset = round(peak_idx / 2);
    x = max_offset:size(ref_wf, 1) + max_offset - 1;
    buffer_wf = size(ref_wf, 1) + 2 * max_offset;
    xq = linspace(1, buffer_wf, buffer_wf * interp_factor);

    interp_wf = interp1(x, ref_wf, xq, "makima");
    interp_wf = interp_wf((max_offset * interp_factor) + 1:((buffer_wf - max_offset - 1) * interp_factor), :);
    col_max = max(abs(interp_wf));
    col_max(col_max == 0) = 1;
    interp_wf = interp_wf ./ col_max;
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
    aligned_wf = aligned_wf((5 * interp_factor):(end - 5 * interp_factor), :);
end
