function [frequency, magnitude, fit_coeff] = computeRegularity(binned_signal, binning, n_peaks)
% COMPUTEREGULARITY  Compute frequency, magnitude and exponential-decay fit of a binned signal.
%
% Computes the dominant frequency and magnitude via FFT, then fits an
% exponential decay to the log-power of the top N peaks. Used identically
% for both single-unit (Unit.getRegularity) and network-level
% (MEArecording.getRegularity) regularity analysis.
%
% INPUTS:
%   binned_signal - (1 x T) binned and normalised activity trace
%   binning       - Bin width in seconds (determines frequency axis)
%   n_peaks       - Number of dominant peaks to use for the decay fit
%
% OUTPUTS:
%   frequency     - Dominant frequency [Hz]
%   magnitude     - FFT magnitude at the dominant frequency
%   fit_coeff     - Exponential decay coefficient (b in fit(x, y, 'exp1'));
%                   NaN if the fit fails

NFFT = length(binned_signal);
F = (0 : 1/NFFT : 1/2-1/NFFT) * (1/binning);
% Apply Hanning window to reduce spectral leakage
win = hanning(NFFT)';
binned_signal = binned_signal .* win;
TEMP = fft(binned_signal, NFFT);
TEMP(1) = 0;
freq_domain = abs(TEMP(1:floor(NFFT/2)));

[magnitude, freq_idx] = max(freq_domain);
frequency = F(freq_idx);

norm_freq_domain = freq_domain / max(freq_domain);
[all_pks, all_locs, ~, all_prom] = findpeaks(norm_freq_domain);
if length(all_pks) >= n_peaks
    [~, prom_order] = sort(all_prom, 'descend');
    sel = prom_order(1:n_peaks);
    [abs_locs, sort_idx] = sort(all_locs(sel));
    p = all_pks(sel(sort_idx));
    abs_locs = abs_locs(:);
    p = p(:);
else
    abs_locs = all_locs(:);
    p = all_pks(:);
end

if length(p) >= 2
    p(p == 0) = eps;  % guard against log10(0) = -Inf from zero-power bins
    log_p = log10(p) - min(log10(p)); % shift to positive for fitting
    try
        f = fit(abs_locs, log_p, 'exp1');
        fit_coeff = f.b;
    catch
        fit_coeff = NaN;
    end
else
    fit_coeff = NaN;
end

end
