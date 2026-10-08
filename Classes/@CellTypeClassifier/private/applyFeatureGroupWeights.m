function X = applyFeatureGroupWeights(X, feat_names, p_umap)
% APPLYFEATUREGROUPWEIGHTS  Dimension-normalised ACG / waveform group scaling.
%
% Columns of each group are multiplied by sqrt(w / n_group), so the group's total
% squared contribution to Euclidean distances is w regardless of how many bins it has
% (1.0 = the whole ACG counts like one scalar feature). Without it, the ~100 ACG bins
% dominate every distance.
%
% INPUTS:
%   X          - (N x F) normalised feature matrix
%   feat_names - (1 x F) feature names matching the columns of X
%   p_umap     - Parameters.UMAP (ACGWeight, WaveformWeight)
    feat_names = string(feat_names);
    assert(numel(feat_names) == size(X, 2), 'CellTypeClassifier:featureNames', ...
        'Got %d feature names for %d columns.', numel(feat_names), size(X, 2));
    is_acg = startsWith(feat_names, "FullACG") | startsWith(feat_names, "ACG");
    is_wf  = startsWith(feat_names, "Waveform") | startsWith(feat_names, "ReferenceWaveform");
    if any(is_acg)
        X(:, is_acg) = X(:, is_acg) * sqrt(p_umap.ACGWeight / sum(is_acg));
    end
    if any(is_wf)
        X(:, is_wf) = X(:, is_wf) * sqrt(p_umap.WaveformWeight / sum(is_wf));
    end
end
