function [W, info] = buildWNNGraph(X, is_acg, k, k_modality)
% BUILDWNNGRAPH  Weighted-nearest-neighbour graph of two feature modalities (ACG, waveform).
%
%   [W, info] = buildWNNGraph(X, is_acg, k)
%
% Multimodal graph as in Seurat's WNN (Hao et al. 2021, Cell) and PhysMAP (Lee et al.
% 2026, Nat Commun): every unit gets its own modality weights, from how well its
% neighbours in each modality predict its features in that modality compared with its
% neighbours in the other modality. Units whose ACG is noisy (few spikes) or whose
% waveform is uninformative then lean on the other modality, and the fixed ACG:waveform
% balance (UMAP.ACGWeight / WaveformWeight) is no longer needed: each modality's
% distances are scaled by their own per-unit bandwidth.
%
% INPUTS:
%   X          - (N x F) normalised features
%   is_acg     - (1 x F) logical, true for ACG columns (the rest are waveform columns)
%   k          - neighbours per unit in the final graph
%   k_modality - neighbours for the modality weights and bandwidths (default 20, as Seurat)
%
% OUTPUTS:
%   W    - (N x N) sparse symmetric graph (fuzzy union of the directed kNN affinities,
%          as in buildKNNGraph)
%   info - struct: weights (N x 2, [ACG waveform]) and the modality names

    if nargin < 4 || isempty(k_modality), k_modality = 20; end
    N = size(X, 1);
    k = min(k, N - 1);
    k_modality = min(k_modality, N - 1);
    % Each modality as in Seurat: a PCA reduction (up to 30 components) with every unit's
    % vector L2-normalised, so distances are bounded and comparable between modalities.
    % On the raw high-dimensional features the cross-modal prediction is almost always
    % ~0, the within / cross ratios explode and the softmax degenerates to choosing one
    % modality per unit.
    mods = {embedModality(X(:, is_acg)), embedModality(X(:, ~is_acg))};
    M = numel(mods);

    nn = cell(1, M); d1 = cell(1, M); sigma = cell(1, M);
    for m = 1:M
        [idx, d] = knnsearch(mods{m}, mods{m}, 'K', max(k, k_modality) + 1);
        idx = idx(:, 2:end); d = d(:, 2:end);              % drop the unit itself
        nn{m} = idx;
        d1{m} = d(:, 1);                                    % distance to the nearest neighbour
        sigma{m} = max(d(:, k_modality) - d1{m}, 1e-12);    % bandwidth (Seurat: k-th minus 1st NN)
    end
    kern = @(dist, m, i) exp(-max(dist - d1{m}(i), 0) ./ sigma{m}(i));

    % Modality weights: within-modality vs cross-modality prediction of each unit
    ratio = zeros(N, M);
    for m = 1:M
        other = 3 - m;
        within = zeros(N, size(mods{m}, 2)); cross = within;
        for i = 1:N
            within(i, :) = mean(mods{m}(nn{m}(i, 1:k_modality), :), 1);
            cross(i, :)  = mean(mods{m}(nn{other}(i, 1:k_modality), :), 1);
        end
        all_i = (1:N)';
        th_w = kern(vecnorm(mods{m} - within, 2, 2), m, all_i);
        th_c = kern(vecnorm(mods{m} - cross, 2, 2), m, all_i);
        ratio(:, m) = th_w ./ (th_c + 1e-4);
    end
    % Linear normalisation of the ratios, as PhysMAP. Seurat's softmax suits strongly
    % coupled modalities (RNA, protein); ACG and waveform are only weakly coupled, the
    % ratios are large for both and the softmax saturates at 0 / 1.
    weights = ratio ./ max(sum(ratio, 2), 1e-12);

    % Weighted affinities to the union of both modalities' neighbours; keep the top k
    rows = zeros(N * k, 1); cols = rows; vals = rows; p = 0;
    for i = 1:N
        cand = unique([nn{1}(i, 1:k), nn{2}(i, 1:k)]);
        s = zeros(1, numel(cand));
        for m = 1:M
            dist = vecnorm(mods{m}(cand, :) - mods{m}(i, :), 2, 2)';
            s = s + weights(i, m) * kern(dist, m, i);
        end
        [s_top, o] = maxk(s, k);
        n_top = numel(o);
        rows(p + (1:n_top)) = i; cols(p + (1:n_top)) = cand(o); vals(p + (1:n_top)) = s_top;
        p = p + n_top;
    end
    Wd = sparse(rows(1:p), cols(1:p), vals(1:p), N, N);
    W  = Wd + Wd' - Wd .* Wd';                              % fuzzy union, as buildKNNGraph
    W(1:N+1:end) = 0;
    info = struct('weights', weights, 'modalities', ["ACG", "Waveform"]);
end

function E = embedModality(Xm, n_pc)
% PCA reduction of one modality (up to n_pc components), rows L2-normalised.
    if nargin < 2, n_pc = 30; end
    Xm(~isfinite(Xm)) = 0;
    Xm = Xm - mean(Xm, 1);
    n_pc = min([n_pc, size(Xm, 2), size(Xm, 1) - 1]);
    [U, S, ~] = svd(Xm, 'econ');
    E = U(:, 1:n_pc) * S(1:n_pc, 1:n_pc);
    E = E ./ max(vecnorm(E, 2, 2), 1e-12);
end
