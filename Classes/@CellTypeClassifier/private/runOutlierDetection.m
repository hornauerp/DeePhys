function [tf_outlier, info] = runOutlierDetection(candidate_coords, contamination, gmm_sep)
% RUNOUTLIERDETECTION  Isolation forest with fixed or auto contamination.
%   Fixed mode (contamination is numeric): iforest with that ContaminationFraction.
%   Auto mode (contamination == "auto"): iforest without ContaminationFraction,
%     then fit 2-component GMM on anomaly scores. Flag the high-score component
%     only if the two components are well-separated (> gmm_sep pooled stds apart).

    n = size(candidate_coords, 1);

    if isstring(contamination) || ischar(contamination)
        is_auto = lower(string(contamination)) == "auto";
    else
        is_auto = false;
    end

    if ~is_auto
        [~, tf_outlier] = iforest(candidate_coords, ...
            'ContaminationFraction', contamination);
        info.mode   = sprintf('fixed=%.2f', contamination);
        info.detail = '';
        return
    end

    if n < 10
        tf_outlier = false(n, 1);
        info.mode   = 'auto';
        info.detail = ' — skipped (n < 10)';
        return
    end

    [forest, ~] = iforest(candidate_coords);
    [~, scores] = isanomaly(forest, candidate_coords);

    try
        gmm = fitgmdist(scores(:), 2, 'RegularizationValue', 1e-6, ...
            'Options', statset('MaxIter', 200, 'TolFun', 1e-6));
    catch
        tf_outlier = false(n, 1);
        info.mode   = 'auto';
        info.detail = ' — GMM fit failed, no outliers removed';
        return
    end

    [~, outlier_comp] = max(gmm.mu);

    pooled_std = sqrt(mean(gmm.Sigma(:)));
    if pooled_std == 0; pooled_std = 1; end
    separation = abs(gmm.mu(1) - gmm.mu(2)) / pooled_std;

    if separation >= gmm_sep
        posterior  = gmm.posterior(scores(:));
        tf_outlier = posterior(:, outlier_comp) > 0.5;
        info.mode   = 'auto';
        info.detail = sprintf(' — GMM separation=%.1f', separation);
    else
        tf_outlier = false(n, 1);
        info.mode   = 'auto';
        info.detail = sprintf(' — no bimodal separation (%.1f < %.1f)', ...
            separation, gmm_sep);
    end
end
