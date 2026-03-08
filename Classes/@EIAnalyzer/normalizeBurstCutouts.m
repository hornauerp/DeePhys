function [norm_bursts, norm_ie] = normalizeBurstCutouts(eia)
% NORMALIZEBURSTCUTOUTS  Normalise burst cutouts across cultures for cross-culture comparison.
%
% Computes the mean burst profile per culture, normalises within each culture,
% and returns matrices suitable for cross-culture visualisation (e.g. heatmaps).
% Only the second half of the burst window (peak and decay) is returned to avoid
% edge artefacts from cross-correlation alignment.
%
% Requires eia.BurstCutouts (run extractBurstCutouts first).
%
% OUTPUTS:
%   norm_bursts - (N_cultures x N_bins) normalised mean total burst activity
%   norm_ie     - (N_cultures x N_bins) normalised mean I/E ratio during bursts

arguments
    eia EIAnalyzer
end

assert(~isempty(eia.BurstCutouts), 'Run extractBurstCutouts() before normalizeBurstCutouts()');

p           = eia.Parameters.BurstDetection;
burst_win   = (p.PeakCutout + 1) : (2*p.PeakCutout + 1);  % peak to end
n_cultures  = length(eia.BurstCutouts);
n_bins      = length(burst_win);

norm_bursts = nan(n_cultures, n_bins);
norm_ie     = nan(n_cultures, n_bins);

for c = 1:n_cultures
    bc = eia.BurstCutouts(c);
    if isempty(bc.total); continue; end

    mean_total = mean(bc.total, 1, 'omitnan');
    mean_i     = mean(bc.i_mat, 1, 'omitnan');

    segment       = mean_total(burst_win);
    norm_bursts(c,:) = segment ./ max(segment);

    ie_segment    = mean_i(burst_win) ./ mean_total(burst_win);
    norm_ie(c,:)  = ie_segment ./ max(ie_segment);
end
end
