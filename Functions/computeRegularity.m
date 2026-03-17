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
abs_locs = [];  % Absolute bin positions in the original freq_domain
p = [];
offset = 0;  % Track cumulative offset from truncation
for i = 1:n_peaks
    if length(norm_freq_domain) < 10
        break
    end
    [pks, locs, ~, ~] = findpeaks(norm_freq_domain, 'NPeaks', 1, 'SortStr', 'descend');
    if isempty(pks)
        break
    end
    abs_locs = [abs_locs; offset + locs]; %#ok<AGROW>  % Convert to absolute position
    p = [p; pks]; %#ok<AGROW>
    offset = offset + locs - 1;  % Truncation removes indices 1..locs-1
    norm_freq_domain = norm_freq_domain(locs:end);
end

if length(p) >= 2
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
