function [t, info] = inferExternalCutoff(p, min_dbic)
% INFEREXTERNALCUTOFF  Label-free decision threshold on P(inhibitory) of an external dataset.
%
% Within the training domain, 0.5 is the calibrated operating point (balanced seeds). An
% external dataset that attaches to the training graph through a biased set of units (e.g.
% units of another preparation that attach mostly through one class) keeps its ranking but
% shifts all scores, so 0.5 misclassifies one class wholesale. Here the threshold is taken
% from the external scores alone: a two-component Gaussian mixture on logit(p);
% the threshold is where both components are equally probable. Used only when the
% distribution is clearly bimodal (BIC of 1 vs 2 components differs by >= min_dbic);
% otherwise 0.5 is kept.
%
% INPUTS:
%   p        - P(inhibitory) of the external units (NaN ignored)
%   min_dbic - minimum BIC improvement of 2 over 1 component (default 10)
% OUTPUTS:
%   t    - threshold (units with p > t are inhibitory)
%   info - struct: method, dBIC, inhibitory_share (weight of the high component), bimodal
    if nargin < 2, min_dbic = 10; end
    p = p(:); p = p(~isnan(p));
    info = struct('method', "gmm", 'dBIC', NaN, 'inhibitory_share', NaN, 'bimodal', false);
    t = 0.5;
    if numel(p) < 20
        info.method = "fixed (too few units)";
        return
    end
    pc = min(max(p, 1e-4), 1 - 1e-4);
    x = log(pc ./ (1 - pc));
    opts = statset('MaxIter', 1000);
    rs = rng; rng(0, 'twister'); restore = onCleanup(@() rng(rs));
    g1 = fitgmdist(x, 1, 'Options', opts);
    g2 = fitgmdist(x, 2, 'Replicates', 20, 'RegularizationValue', 1e-4, 'Options', opts);
    info.dBIC = g1.BIC - g2.BIC;
    [~, hi] = max(g2.mu);
    info.inhibitory_share = g2.ComponentProportion(hi);
    if info.dBIC < min_dbic
        info.method = "fixed (scores not bimodal)";
        return
    end
    grid = linspace(min(g2.mu), max(g2.mu), 2000)';
    post = posterior(g2, grid);
    [~, k] = min(abs(post(:, hi) - 0.5));
    t = 1 / (1 + exp(-grid(k)));
    info.bimodal = true;
end
