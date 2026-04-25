function coords = prepareIforestInput(X_feat, local_indices, ...
        reduction, subset_global_idx, domain, p_outlr)
% PREPAREIFORESTINPUT  Get coordinates for iforest based on domain.
%   "umap" — index into the full UMAP embedding using global indices
%   "pca"  — PCA on candidate features (default)

    if domain == "umap"
        global_idx = subset_global_idx(local_indices);
        coords = reduction(global_idx, :);
    else
        candidate_features = X_feat(local_indices, :);
        col_var = var(candidate_features, 0, 1);
        candidate_features_clean = candidate_features(:, col_var > 1e-10);
        n_pca = min([p_outlr.NPCAComponents, ...
                     size(candidate_features_clean, 1) - 1, ...
                     size(candidate_features_clean, 2)]);
        if n_pca >= 1
            [~, coords] = pca(candidate_features_clean, ...
                'NumComponents', n_pca, 'Algorithm', 'eig');
        else
            coords = candidate_features_clean;
        end
    end
end
