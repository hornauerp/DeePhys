function [clean_global, outlier_mask] = graphPurityOutlier( ...
        reduction, candidate_global, is_responsive_global, k, sigma_thresh)
% GRAPHPURITYOUTLIER  Remove responsive units whose UMAP neighbourhood is
%   dominated by non-responsive units (likely mislabelled or boundary units).
%
% For each unit in candidate_global, computes the fraction of its k nearest
% neighbours (in the full UMAP embedding) that are also responsive. Units
% whose purity is more than sigma_thresh robust-z-scores below the median
% purity are flagged as outliers.

    if isempty(candidate_global)
        clean_global = candidate_global;
        outlier_mask = false(1, 0);
        return
    end

    N_all = size(reduction, 1);
    k_act = min(k, N_all - 1);
    M     = numel(candidate_global);

    if k_act < 1 || M == 0
        clean_global = candidate_global;
        outlier_mask = false(1, M);
        return
    end

    [nn_all, ~] = knnsearch(reduction, reduction(candidate_global, :), 'K', k_act + 1);
    nn_all = nn_all(:, 2:end);   % exclude self

    purity = mean(is_responsive_global(nn_all), 2);

    med_p = median(purity);
    mad_p = median(abs(purity - med_p));
    if mad_p == 0
        outlier_mask = false(1, M);
        clean_global = candidate_global;
        return
    end

    robust_z  = (purity - med_p) / (1.4826 * mad_p);
    outlier_mask = (robust_z < -sigma_thresh)';
    clean_global = candidate_global(~outlier_mask);
end
