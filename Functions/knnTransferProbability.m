function p = knnTransferProbability(X_ref, p_ref, X_query, k)
% KNNTRANSFERPROBABILITY  P(inhibitory) of query units from their nearest reference units.
%
%   p = knnTransferProbability(X_ref, p_ref, X_query, k)
%
% Each query unit gets the mean p_ref of its k nearest reference units (Euclidean), weighted
% by 1 / distance (distances below 1e-9 count as 1e-9). Reference units with NaN p_ref are
% left out; without any, p is NaN. Used by CellTypeClassifier.classifyExternalUnits
% (Method "knn").
%
% INPUTS:
%   X_ref   - (N_ref x F) reference features
%   p_ref   - (N_ref x 1) reference P(inhibitory)
%   X_query - (N_query x F) query features, normalised like X_ref
%   k       - number of neighbours (capped at the number of reference units with p_ref)
% OUTPUT:
%   p       - (N_query x 1) P(inhibitory)

arguments
    X_ref   (:,:) double
    p_ref   (:,1) double
    X_query (:,:) double
    k       (1,1) double {mustBePositive, mustBeInteger}
end
have = isfinite(p_ref);
p = nan(size(X_query, 1), 1);
if ~any(have) || isempty(X_query), return; end
k = min(k, nnz(have));
[idx, d] = knnsearch(X_ref(have, :), X_query, 'K', k);
w = 1 ./ max(d, 1e-9);
pr = p_ref(have);
p = sum(w .* reshape(pr(idx), size(idx)), 2) ./ sum(w, 2);
end
