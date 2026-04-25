function [X_all, norm_params] = computeGlobalNormalization(nf, p_umap, ph)
% COMPUTEGLOBALNORMALIZATION  Global z-score, NaN removal, scale, optional feature selection.
%
% Fits on ALL unique units (transductive setting). Returns the normalized matrix
% and the NormalizationParams struct consumed by classifyUnits / classifyExternalUnits.

[X_all, mu_global, sigma_global] = normalize(nf.X_pergroup);

nan_cols = any(isnan(X_all), 1);
X_all(:, nan_cols) = [];
scale = max(abs(X_all), [], 1);
scale(scale == 0) = 1;
X_all = X_all ./ scale;

feat_names_trimmed = nf.feat_names(~nan_cols);

if p_umap.FeatureSelection && size(X_all, 2) > 1
    feat_var   = var(X_all, 0, 1);
    var_thresh = prctile(feat_var, p_umap.MinVariancePercentile);
    low_var    = feat_var <= var_thresh;

    % Graph-based correlation removal: iteratively remove the feature with the
    % most high-correlation neighbors (order-independent, avoids systematic bias).
    R = abs(corrcoef(X_all));
    R(logical(eye(size(R)))) = 0;
    high_corr = false(1, size(X_all, 2));
    adj = R > p_umap.MaxCorrelation;
    while any(adj(:))
        n_edges = sum(adj & ~high_corr' & ~high_corr, 2);
        n_edges(high_corr) = 0;
        [~, worst] = max(n_edges);
        high_corr(worst) = true;
        adj(worst, :) = false;
        adj(:, worst) = false;
    end

    remove_feat = low_var | high_corr;
    feature_selection_mask = ~remove_feat;
    X_all = X_all(:, feature_selection_mask);
    feat_names_trimmed = feat_names_trimmed(feature_selection_mask);
    fprintf('Feature selection: removed %d low-var + %d redundant = %d total (kept %d)\n', ...
        sum(low_var), sum(high_corr & ~low_var), sum(remove_feat), size(X_all, 2));
else
    feature_selection_mask = true(1, numel(feat_names_trimmed));
end

% Compute actual waveform window params so classifyExternalUnits can verify
% that external features were built with the same window dimensions.
n_wf_features    = sum(startsWith(nf.feat_names, "Waveform"));
target_sr        = ph.WaveformTargetSamplingRate;
dt_ms            = 1000 / target_sr;
n_pre_requested  = round(ph.WaveformPreTrough  * target_sr / 1000);
n_post_requested = round(ph.WaveformPostTrough * target_sr / 1000);
n_requested      = n_pre_requested + n_post_requested + 1;

if n_wf_features == n_requested
    wf_pre_ms  = ph.WaveformPreTrough;
    wf_post_ms = ph.WaveformPostTrough;
else
    % Trim mode narrowed the window; back-compute actual pre/post from features.
    ratio             = n_pre_requested / (n_pre_requested + n_post_requested);
    n_out_pre_actual  = round((n_wf_features - 1) * ratio);
    n_out_post_actual = n_wf_features - 1 - n_out_pre_actual;
    wf_pre_ms         = n_out_pre_actual  * dt_ms;
    wf_post_ms        = n_out_post_actual * dt_ms;
end

norm_params = struct( ...
    'mu_global',              mu_global, ...
    'sigma_global',           sigma_global, ...
    'nan_cols',               nan_cols, ...
    'scale',                  scale, ...
    'feature_selection_mask', feature_selection_mask, ...
    'feat_names_trimmed',     feat_names_trimmed, ...
    'wf_pre_trough_ms',       wf_pre_ms, ...
    'wf_post_trough_ms',      wf_post_ms);
end
