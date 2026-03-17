function Burst = detectBurstsLogISI(obj, min_spikes)
% DETECTBURSTLOGISI  Detect bursts using the logISI method (Pasquale et al., 2010).
%
%   Fits a two-component Gaussian mixture model to the log10(ISI) distribution
%   to find a data-driven threshold separating intra-burst from inter-burst ISIs.
%   This avoids the fixed-threshold assumption of the ISIN method.
%
% INPUTS:
%   obj        - MEArecording with populated Spikes
%   min_spikes - Minimum number of spikes per burst (default: 3)
%
% OUTPUTS:
%   Burst      - struct with T_start and T_end vectors
%
% REFERENCE:
%   Pasquale V, Martinoia S, Bhalla US (2010). A self-adaptive approach
%   to burst detection. J Neurosci Methods, 191(1):88-96.

    arguments
        obj MEArecording
        min_spikes (1,1) double = 3
    end

    spike_times = obj.Spikes.Times;

    if length(spike_times) < min_spikes * 2
        Burst.T_start = [];
        Burst.T_end = [];
        return
    end

    isi = diff(spike_times);
    log_isi = log10(isi);

    % Remove extreme outliers before fitting
    log_isi = log_isi(isfinite(log_isi));
    if length(log_isi) < 20
        Burst.T_start = [];
        Burst.T_end = [];
        return
    end

    % Fit 2-component Gaussian mixture to log(ISI) distribution
    try
        gm = fitgmdist(log_isi(:), 2, 'RegularizationValue', 1e-5, ...
            'Options', statset('MaxIter', 500, 'TolFun', 1e-6), 'Replicates', 3);
    catch
        % If GMM fails, fall back to void result
        Burst.T_start = [];
        Burst.T_end = [];
        return
    end

    % Identify which component has the smaller mean (intra-burst = short ISIs)
    [~, sort_idx] = sort(gm.mu);
    intra_comp = sort_idx(1);  % shorter ISIs = intra-burst
    inter_comp = sort_idx(2);  % longer ISIs = inter-burst

    % Check bimodality: if means are not well-separated, no clear burst structure
    separation = abs(diff(gm.mu)) / mean(sqrt(gm.Sigma(:)));
    if separation < 1.0
        % Unimodal — no clear burst/non-burst separation
        Burst.T_start = [];
        Burst.T_end = [];
        return
    end

    % Find threshold: intersection point of the two Gaussians
    % Search between the two means
    mu1 = gm.mu(intra_comp); sigma1 = sqrt(gm.Sigma(intra_comp));
    mu2 = gm.mu(inter_comp); sigma2 = sqrt(gm.Sigma(inter_comp));
    w1 = gm.ComponentProportion(intra_comp);
    w2 = gm.ComponentProportion(inter_comp);

    % Fine grid search for the crossing point
    x_range = linspace(mu1, mu2, 1000);
    pdf1 = w1 * normpdf(x_range, mu1, sigma1);
    pdf2 = w2 * normpdf(x_range, mu2, sigma2);
    [~, cross_idx] = min(abs(pdf1 - pdf2));
    log_threshold = x_range(cross_idx);
    threshold = 10^log_threshold;  % Convert back from log10 to seconds

    % Classify ISIs as intra-burst (below threshold) or inter-burst
    is_intra = isi <= threshold;

    % Build bursts from consecutive intra-burst ISIs
    burst_starts = [];
    burst_ends = [];
    in_burst = false;
    b_start = 0;

    for i = 1:length(is_intra)
        if is_intra(i) && ~in_burst
            b_start = i;
            in_burst = true;
        elseif ~is_intra(i) && in_burst
            % End of burst: ISI i is inter-burst, so last intra-burst spike is i
            % (ISI b_start..i-1 are intra-burst, spanning spikes b_start..i)
            n_spikes_in_burst = i - b_start + 1;
            if n_spikes_in_burst >= min_spikes
                burst_starts = [burst_starts; spike_times(b_start)]; %#ok<AGROW>
                burst_ends = [burst_ends; spike_times(i)]; %#ok<AGROW> % spike_times(i) is the last spike connected by intra-burst ISIs
            end
            in_burst = false;
        end
    end

    % Handle burst extending to end of recording
    if in_burst
        n_spikes_in_burst = length(is_intra) - b_start + 2;
        if n_spikes_in_burst >= min_spikes
            burst_starts = [burst_starts; spike_times(b_start)];
            burst_ends = [burst_ends; spike_times(end)];
        end
    end

    Burst.T_start = burst_starts;
    Burst.T_end = burst_ends;
    Burst.LogISI_Threshold = threshold;
    Burst.GMModel = gm;
end
