function [T, post, p, mdl_all] = estimateClassComposition(Z, train_rows, y, groups, opts)
% ESTIMATECLASSCOMPOSITION  Per-group class fraction and per-unit identity from calibrated
% probabilities (quantification under prior shift).
%
%   [T, post, p] = estimateClassComposition(Z, train_rows, y, groups)
%
% 1. Ridge logistic regression on the training units (labels y: 2 = positive class, e.g.
%    inhibitory). Training units get out-of-fold probabilities (cross-fitting with folds by
%    group, i.e. culture), all others the model fitted on every training unit.
% 2. Per group (culture / well): EM re-estimation of the group's positive fraction and of the
%    unit posteriors under that fraction (adjustClassPrior; Saerens et al. 2002), with a
%    bootstrap interval over the group's units.
% Unit identity is kept: post is each unit's probability given its features and its group's
% estimated composition; labels follow from post > 0.5.
%
% INPUTS:
%   Z          - (N x D) features (e.g. principal components of the harmonised features)
%   train_rows - indices (or logical) of the training units
%   y          - labels of the training units (1 / 2), same order as train_rows
%   groups     - (N x 1) group of each unit
%   opts.NFolds (5), opts.NBoot (200), opts.Lambda ([] = 1 / n_train), opts.Seed (1)
% OUTPUTS:
%   T    - table: Group, NUnits, Fraction, CILow, CIHigh, Counted (share p > 0.5),
%          MeanProbability
%   post - (N x 1) posterior after the group's prior adjustment
%   p    - (N x 1) classifier probability at the training prior
%   mdl_all - the logistic model fitted on all training units (scores new rows of Z)
arguments
    Z          double
    train_rows
    y
    groups
    opts.NFolds (1,1) double = 5
    opts.NBoot  (1,1) double = 200
    opts.Lambda double = []
    opts.Seed   (1,1) double = 1
end
N = size(Z, 1);
if islogical(train_rows), train_rows = find(train_rows); end
train_rows = train_rows(:); y = double(y(:)) == 2;
assert(numel(train_rows) == numel(y), 'estimateClassComposition:labels', ...
    'train_rows and y must have the same length.');
assert(any(y) && any(~y), 'estimateClassComposition:oneClass', 'Training labels need both classes.');
groups = string(groups(:));
lambda = opts.Lambda;
if isempty(lambda), lambda = 1 / numel(train_rows); end
rng(opts.Seed, 'twister');

fit = @(rows, yy) fitclinear(Z(rows, :), yy, 'Learner', 'logistic', ...
    'Regularization', 'ridge', 'Lambda', lambda);
prob = @(mdl, rows) localProb(mdl, Z(rows, :));

% -- probabilities: cross-fitted for training units, full model for the rest -----------
p = nan(N, 1);
mdl_all = fit(train_rows, y);
p(:) = prob(mdl_all, (1:N)');
% Folds by group: units of one culture are correlated, so a unit-level split would score each
% training unit with a model that saw its culture mates (optimistic probabilities).
[~, ~, gi] = unique(groups(train_rows));
fold = zeros(numel(train_rows), 1);
if max(gi) >= 2
    k = min(opts.NFolds, max(gi));
    fold_of_group = mod(randperm(max(gi)), k) + 1;
    fold = fold_of_group(gi);
    fold = fold(:);
else
    k = min(opts.NFolds, min(sum(y), sum(~y)));   % a single group: unit-level folds
    if k >= 2
        cv = cvpartition(y, 'KFold', k);
        for f = 1:k, fold(test(cv, f)) = f; end
    end
end
for f = 1:max(fold)
    te = fold == f; tr = ~te;
    if ~any(y(tr)) || all(y(tr)), continue, end   % both classes needed; keep the full-model p
    mdl = fit(train_rows(tr), y(tr));
    p(train_rows(te)) = prob(mdl, train_rows(te));
end
pi_tr = mean(y);

% -- per group: EM prior adjustment with a bootstrap interval ------------------------
[G, gname] = findgroups(groups);
ng = numel(gname);
frac = nan(ng, 1); lo = frac; hi = frac; counted = frac; meanp = frac; n_units = zeros(ng, 1);
post = nan(N, 1);
for g = 1:ng
    rows = find(G == g);
    pg = p(rows);
    [frac(g), post(rows)] = adjustClassPrior(pg, pi_tr);
    boot = nan(opts.NBoot, 1);
    for b = 1:opts.NBoot
        boot(b) = adjustClassPrior(pg(randi(numel(pg), numel(pg), 1)), pi_tr);
    end
    lo(g) = prctile(boot, 2.5); hi(g) = prctile(boot, 97.5);
    counted(g) = mean(pg > 0.5); meanp(g) = mean(pg); n_units(g) = numel(rows);
end
T = table(gname, n_units, frac, lo, hi, counted, meanp, 'VariableNames', ...
    {'Group', 'NUnits', 'Fraction', 'CILow', 'CIHigh', 'Counted', 'MeanProbability'});
end

function q = localProb(mdl, X)
    [~, s] = predict(mdl, X);
    q = s(:, 2);
end
