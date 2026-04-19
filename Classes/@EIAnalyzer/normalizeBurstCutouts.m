function normalizeBurstCutouts(eia)
% NORMALIZEBURSTCUTOUTS  Normalise burst cutouts across cultures for cross-culture comparison.
%
% Computes the mean burst profile per culture from BurstCutouts.accepted,
% normalises within each culture, and assembles a N_cultures x N_bins matrix.
%
% NormalizationWindow controls which portion of the burst window is stored:
%   "full"        — entire burst window including pre-peak ramp-up (default)
%   "peak_to_end" — peak and decay only (bins PeakCutout+1 onward)
%
% The full window is preferred because the pre-peak ramp-up carries information
% about burst initiation dynamics. NaN-padded edges from xcorr alignment are
% handled by omitnan means, so alignment artefacts do not distort the result.
%
% Requires eia.BurstCutouts (run extractBurstCutouts first).
% Sets: eia.NormalizedCutouts with fields:
%   .bursts     — (N_cultures x N_bins) normalized mean total burst profile
%   .inh_frac   — (N_cultures x N_bins) normalized inhibitory fraction (I/total),
%                 i.e. the proportion of total activity that is inhibitory during bursts.
%                 This is NOT the I/E ratio; it is normalized to its own maximum within
%                 each culture so only burst shape is comparable across cultures.
%   .burst_counts — (N_cultures x 1) number of accepted bursts per culture

arguments
    eia EIAnalyzer
end

assert(~isempty(eia.BurstCutouts), 'Run extractBurstCutouts() before normalizeBurstCutouts()');

p          = eia.Parameters.BurstDetection;
n_cultures = numel(eia.BurstCutouts);

switch p.NormalizationWindow
    case "peak_to_end"
        burst_win = (p.PeakCutout + 1) : (2*p.PeakCutout + 1);
    case "full"
        burst_win = 1 : (2*p.PeakCutout + 1);
    otherwise
        warning('EIAnalyzer:normalizeBurstCutouts', ...
            'Unknown NormalizationWindow "%s" — using "full".', p.NormalizationWindow);
        burst_win = 1 : (2*p.PeakCutout + 1);
end

n_bins       = numel(burst_win);
norm_bursts  = nan(n_cultures, n_bins);
norm_inh_frac = nan(n_cultures, n_bins);
burst_counts = zeros(n_cultures, 1);

for c = 1:n_cultures
    fprintf('  normalizeBurstCutouts: culture %d/%d\n', c, n_cultures);
    bc = eia.BurstCutouts(c).accepted;
    if isempty(bc.total); continue; end
    burst_counts(c) = size(bc.total, 1);

    mean_total = mean(bc.total, 1, 'omitnan');
    mean_i     = mean(bc.i_mat, 1, 'omitnan');

    segment = mean_total(burst_win);
    seg_max = max(segment);
    if seg_max == 0; continue; end

    % Inhibitory fraction: proportion of total activity that is inhibitory (I/total).
    % Note: this is NOT the I/E ratio. Normalized to its own maximum within each culture.
    inh_frac = mean_i(burst_win) ./ mean_total(burst_win);
    inh_frac(~isfinite(inh_frac)) = 0;
    inh_frac_max = max(inh_frac);
    if inh_frac_max == 0; continue; end

    norm_bursts(c,:)   = segment  ./ seg_max;
    norm_inh_frac(c,:) = inh_frac ./ inh_frac_max;
end

eia.NormalizedCutouts = struct('bursts', norm_bursts, 'inh_frac', norm_inh_frac, 'burst_counts', burst_counts);
fprintf('Normalised burst cutouts for %d cultures\n', n_cultures);
end
