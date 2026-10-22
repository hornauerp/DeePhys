function [X, fallback, stats] = zscoreByGroup(X, groups, is_ref, min_ref, opts)
% ZSCOREBYGROUP  Z-score the rows of X within each group, optionally with reference statistics.
%
%   X = zscoreByGroup(X, groups)
%   [X, fallback, stats] = zscoreByGroup(X, groups, is_ref, min_ref)
%   X = zscoreByGroup(X_new, groups_new, Stats=stats)     % apply stored statistics
%
% Each group is centred and scaled with the mean / SD of its reference rows (is_ref) only.
% With bridge rows of equal composition in every group (e.g. the same culture condition on
% every plate), technical offsets between groups are removed while differences in the
% composition of the other rows are kept; plain per-group z-scoring would remove both.
%
% INPUTS:
%   X       - (N x F) features
%   groups  - (N x 1) group labels (numeric, string, categorical)
%   is_ref  - (N x 1) logical reference rows (default: all rows)
%   min_ref - groups with fewer reference rows use all their rows (default 20)
%   ConstantAs - value for columns without variance in a group: 0 (default) or NaN (as
%                MATLAB's normalize, which the CellTypeClassifier relies on: NaN columns
%                are removed afterwards)
%   Stats   - statistics returned by an earlier call: new rows (e.g. units with resampled
%             spike trains) are standardised with the mean / SD of their group instead of
%             their own; rows of groups not in Stats become NaN
% OUTPUTS:
%   X        - standardised features
%   fallback - group indices (order of first appearance) that fell back to all rows
%   stats    - struct: Group (G x 1 string), Mu, SD (G x F; SD 0 = constant column)
arguments
    X       double
    groups
    is_ref  logical = true(size(X, 1), 1)
    min_ref (1,1) double = 20
    opts.ConstantAs (1,1) double = 0
    opts.Stats struct = struct([])
end
[gname, ~, iG] = unique(string(groups(:)), 'stable');
fallback = zeros(1, 0);
if ~isempty(opts.Stats)
    [hit, loc] = ismember(gname, opts.Stats.Group);
    for g = 1:numel(gname)
        mask = iG == g;
        if ~hit(g), X(mask, :) = NaN; continue; end
        mu = opts.Stats.Mu(loc(g), :); sd = opts.Stats.SD(loc(g), :);
        flat = sd == 0; sd(flat) = 1;
        X(mask, :) = (X(mask, :) - mu) ./ sd;
        X(mask, flat) = opts.ConstantAs;
    end
    stats = opts.Stats;
    return
end
is_ref = is_ref(:);
stats = struct('Group', gname, 'Mu', nan(numel(gname), size(X, 2)), 'SD', nan(numel(gname), size(X, 2)));
for g = 1:max(iG)
    mask  = iG == g;
    rmask = mask & is_ref;
    if sum(rmask) < min_ref
        if ~all(is_ref), fallback(end+1) = g; end %#ok<AGROW>
        rmask = mask;
    end
    mu = mean(X(rmask, :), 1, "omitnan");     % missing rows stay missing
    sd = std(X(rmask, :), 0, 1, "omitnan");
    stats.Mu(g, :) = mu; stats.SD(g, :) = sd;
    flat = sd == 0;
    sd(flat) = 1;
    X(mask, :) = (X(mask, :) - mu) ./ sd;
    X(mask, flat) = opts.ConstantAs;
end
end
