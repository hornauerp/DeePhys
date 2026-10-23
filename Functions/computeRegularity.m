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
%   magnitude     - FFT magnitude at the dominant frequency over the mean magnitude of the
%                   spectrum. Free of the raw magnitude's growth with duration, but as the
%                   maximum over all frequency bins it still rises slowly with their number:
%                   aperiodic (Poisson) signals give about 3-4, more in longer recordings.
%                   Compare recordings of equal duration and binning.
%   fit_coeff     - Exponential decay coefficient (b in fit(f, y, 'exp1'), f in Hz);
%                   NaN if the fit fails or the spectral peak is not significant against
%                   bin-shuffled surrogates (no rhythm: the fit would describe noise)

NFFT = length(binned_signal);
F = (0 : 1/NFFT : 1/2-1/NFFT) * (1/binning);
freq_domain = spectrumMagnitude(binned_signal);

[peak_mag, freq_idx] = max(freq_domain);
frequency = F(freq_idx);
% Peak over the mean spectral magnitude: the raw FFT magnitude grows with the number of
% bins, i.e. with recording duration, even for an aperiodic (Poisson) signal.
magnitude = peak_mag / max(mean(freq_domain), eps);

% The decay fit describes the harmonics of a rhythm; for an aperiodic signal it fits noise
% (unbounded coefficients). It is computed only when the spectral peak exceeds that of all
% 19 bin-shuffled surrogates (p <= 0.05; fixed seed, so results are reproducible).
fit_coeff = NaN;
if NFFT < 8 || ~any(binned_signal)
    return
end
rs = RandStream('mt19937ar', 'Seed', 1);
null_mag = zeros(1, 19);
for k = 1:numel(null_mag)
    fd = spectrumMagnitude(binned_signal(randperm(rs, NFFT)));
    null_mag(k) = max(fd) / max(mean(fd), eps);
end
if magnitude <= max(null_mag)
    return
end

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
        % Fit against frequency in Hz (not FFT bin index), so the coefficient does not
        % scale with recording duration.
        f = fit(F(abs_locs)', log_p, 'exp1');
        fit_coeff = f.b;
    catch
        fit_coeff = NaN;
    end
else
    fit_coeff = NaN;
end

end

function freq_domain = spectrumMagnitude(binned_signal)
% Magnitude spectrum without DC. Mean and linear drift are removed before the Hanning window:
% a positive signal times the window has a strong low-frequency lobe, so the "dominant
% frequency" was always the lowest non-DC bin (1/duration).
NFFT = length(binned_signal);
x = detrend(binned_signal(:))' .* hanning(NFFT)';
TEMP = fft(x, NFFT);
TEMP(1) = 0;
freq_domain = abs(TEMP(1:floor(NFFT/2)));
end
