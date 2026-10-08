function [prior, post, n_iter] = adjustClassPrior(p, train_prior, opts)
% ADJUSTCLASSPRIOR  Class fraction of a new group of units under prior shift (EM of
% Saerens, Latinne & Decaestecker 2002, Neural Computation 14:21-41).
%
%   [prior, post] = adjustClassPrior(p, train_prior)
%
% A probabilistic classifier trained at class fraction train_prior returns posteriors p for
% the units of one culture (well). Counting p > 0.5 or averaging p pulls every culture towards
% the training fraction. The EM re-estimates the culture's own fraction q and the posteriors:
%   post_i = q/pi * p_i / (q/pi * p_i + (1-q)/(1-pi) * (1-p_i)),   q = mean(post).
% Requires calibrated posteriors (e.g. logistic regression).
%
% INPUTS:
%   p           - (N x 1) posterior probability of the positive class per unit
%   train_prior - positive-class fraction of the training set
%   opts.Tol, opts.MaxIter, opts.Init (start value of q, default 0.5)
% OUTPUTS:
%   prior  - estimated positive-class fraction of these units
%   post   - (N x 1) adjusted posteriors
%   n_iter - EM iterations used
arguments
    p           (:,1) double {mustBeInRange(p, 0, 1)}
    train_prior (1,1) double {mustBeInRange(train_prior, 0, 1, 'exclusive')}
    opts.Tol     (1,1) double = 1e-6
    opts.MaxIter (1,1) double = 500
    opts.Init    (1,1) double = 0.5
end
q = opts.Init; post = p;
n_iter = 0;
if isempty(p), prior = NaN; return; end
for n_iter = 1:opts.MaxIter
    w1 = q / train_prior * p;
    w0 = (1 - q) / (1 - train_prior) * (1 - p);
    post = w1 ./ max(w1 + w0, realmin);
    q_new = mean(post);
    if abs(q_new - q) < opts.Tol, q = q_new; break; end
    q = q_new;
end
prior = q;
end
