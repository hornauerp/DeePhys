function normalizeBurstCutouts(eia)
% NORMALIZEBURSTCUTOUTS  Normalise burst cutouts across cultures for cross-culture comparison.
%
% Computes the mean burst profile per culture from BurstCutouts.accepted,
% normalises within each culture, and assembles a N_cultures x N_bins matrix.
%
% NormalizationWindow controls which portion of the burst window is stored:
%   "peak_to_end" — second half only (peak and decay); avoids xcorr alignment
%                   artefacts in the pre-peak region (default)
%   "full"        — entire burst window including pre-peak ramp-up
%
% Requires eia.BurstCutouts (run extractBurstCutouts first).
% Sets: eia.NormalizedCutouts  (.bursts and .ie, each N_cultures x N_bins)

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
            'Unknown NormalizationWindow "%s" — using "peak_to_end".', p.NormalizationWindow);
        burst_win = (p.PeakCutout + 1) : (2*p.PeakCutout + 1);
end

n_bins       = numel(burst_win);
norm_bursts  = nan(n_cultures, n_bins);
norm_ie      = nan(n_cultures, n_bins);
burst_counts = zeros(n_cultures, 1);

for c = 1:n_cultures
    bc = eia.BurstCutouts(c).accepted;
    if isempty(bc.total); continue; end
    burst_counts(c) = size(bc.total, 1);

    mean_total = mean(bc.total, 1, 'omitnan');
    mean_i     = mean(bc.i_mat, 1, 'omitnan');

    segment = mean_total(burst_win);
    seg_max = max(segment);
    if seg_max == 0; continue; end

    ie_segment = mean_i(burst_win) ./ mean_total(burst_win);
    ie_segment(~isfinite(ie_segment)) = 0;
    ie_max = max(ie_segment);
    if ie_max == 0; continue; end

    norm_bursts(c,:) = segment    ./ seg_max;
    norm_ie(c,:)     = ie_segment ./ ie_max;
end

eia.NormalizedCutouts = struct('bursts', norm_bursts, 'ie', norm_ie, 'burst_counts', burst_counts);
fprintf('Normalised burst cutouts for %d cultures\n', n_cultures);
end
