function [sigma_ms, info] = selectACGSmoothing(spike_trains, acg_used, bin_size, lag)
% SELECTACGSMOOTHING  Label-free choice of the ACG smoothing width by split-half prediction.
%
% Each spike train is split into alternating 10-s blocks (two halves of the same recording,
% independent noise, same firing state). For each candidate width the smoothed ACG of one
% half predicts the other half (squared error of sum-normalised ACGs, both directions); the
% best width minimises the error over units. The halves have half the spike pairs of the full
% ACG, so the width is scaled by (pairs_half / pairs_used)^(1/5) (kernel-bandwidth rate) to
% the ACGs actually used as features.
%
% INPUTS:
%   spike_trains - cell array of spike-time vectors (s), e.g. the representative units
%   acg_used     - (n_bins x N) ACGs that will be smoothed (for their pair counts), or []
%   bin_size, lag- ACG parameters (s)
% OUTPUTS:
%   sigma_ms - chosen Gaussian width (ms) for acg_used
%   info     - grid, error curve, half-data optimum, scale factor
    grid_ms = [0 0.25 0.5 0.75 1 1.5 2 3 4];
    n = numel(spike_trains);
    use = find(cellfun(@numel, spike_trains) >= 200);
    if numel(use) > 1000
        use = use(round(linspace(1, numel(use), 1000)));
    end
    err = zeros(size(grid_ms)); pairs_half = nan(numel(use), 1);
    nb = round(2 * lag / bin_size) + 1; c = (nb + 1) / 2;
    for i = 1:numel(use)
        st = sort(spike_trains{use(i)}(:));
        blk = mod(floor(st / 10), 2) == 0;
        A = CellTypeClassifier.computeACG(st(blk), bin_size, lag);
        B = CellTypeClassifier.computeACG(st(~blk), bin_size, lag);
        A(c) = 0; B(c) = 0;
        if sum(A) == 0 || sum(B) == 0, continue, end
        pairs_half(i) = (sum(A) + sum(B)) / 2;
        An = A / sum(A); Bn = B / sum(B);
        for g = 1:numel(grid_ms)
            sb = grid_ms(g) / 1000 / bin_size;
            err(g) = err(g) + sum((smoothACG(An, sb) - Bn).^2) + sum((smoothACG(Bn, sb) - An).^2);
        end
    end
    [~, gi] = min(err);
    sigma_half = grid_ms(gi);
    scale = 1;
    if ~isempty(acg_used)
        % pair counts of the SAME units (columns of acg_used match spike_trains)
        pu = acg_used(:, use); pu(c, :) = 0;
        ratio = pairs_half(:) ./ sum(pu, 1)';
        ratio = ratio(isfinite(ratio) & ratio > 0);
        if ~isempty(ratio), scale = median(ratio)^(1/5); end
    end
    sigma_ms = sigma_half * scale;
    info = struct('grid_ms', grid_ms, 'error', err, 'sigma_half_ms', sigma_half, ...
        'scale', scale, 'n_units', sum(~isnan(pairs_half)), 'n_trains', n);
end
