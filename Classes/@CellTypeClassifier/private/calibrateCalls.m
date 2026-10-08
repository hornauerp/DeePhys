function [inc, dec, info] = calibrateCalls(T, is_ctrl, q, alpha, target_fdr)
% CALIBRATECALLS  Responder calls with thresholds set by control cultures.
%
% Increase: LFCc above a threshold and rank-sum p < alpha; decrease: LFCc below a threshold
% and p < alpha. Controls (no DREADD target) define chance; the estimated false-discovery
% rate of a call type is the control call rate divided by the treated call rate.
%
% Thresholds:
%   target_fdr empty/NaN : the q and (1-q) quantiles of control LFCc (fixed tail)
%   target_fdr given     : the least strict threshold whose estimated FDR <= target_fdr
%                          (most responders at a stated purity; strict tails select only the
%                          strongest responders, i.e. a biased subset of the target cells)
%
% INPUTS:
%   T          - table from calibratedResponseTable (LFCc, p)
%   is_ctrl    - (height(T) x 1) logical, row belongs to a control culture
%   q          - control quantile per direction (e.g. 0.01), used without target_fdr
%   alpha      - rank-sum significance level (e.g. 0.01)
%   target_fdr - (optional) target false-discovery rate, e.g. 0.2
% OUTPUTS:
%   inc, dec - (height(T) x 1) logical calls (controls included, for the FDR estimate)
%   info     - thresholds, call rates and estimated FDR
    if nargin < 5, target_fdr = NaN; end
    assert(any(is_ctrl), 'CellTypeClassifier:noControls', ...
        'No control-culture units: control_calibrated responder calls need control cultures (Bootstrap.ControlFilter).');
    ctrl = T.LFCc(is_ctrl);
    tr = ~is_ctrl;
    sig = T.p < alpha;
    if isempty(target_fdr) || isnan(target_fdr)
        thr_inc = quantile(ctrl, 1 - q);
        thr_dec = quantile(ctrl, q);
    else
        thr_inc = fdrThreshold(T.LFCc, sig, is_ctrl, target_fdr, +1);
        thr_dec = -fdrThreshold(-T.LFCc, sig, is_ctrl, target_fdr, +1);
    end
    inc = T.LFCc > thr_inc & sig;
    dec = T.LFCc < thr_dec & sig;
    rate = @(x, m) mean(x(m));
    info = struct('thr_inc', thr_inc, 'thr_dec', thr_dec, 'q', q, 'alpha', alpha, 'target_fdr', target_fdr, ...
        'n_control', sum(is_ctrl), 'n_treated', sum(tr), ...
        'rate_inc_control', rate(inc, is_ctrl), 'rate_inc_treated', rate(inc, tr), ...
        'rate_dec_control', rate(dec, is_ctrl), 'rate_dec_treated', rate(dec, tr));
    info.fdr_inc = min(1, info.rate_inc_control / max(info.rate_inc_treated, eps));
    info.fdr_dec = min(1, info.rate_dec_control / max(info.rate_dec_treated, eps));
end

function thr = fdrThreshold(x, sig, is_ctrl, target, ~)
% Smallest threshold t (scanning candidate values of x among treated significant units)
% such that rate_control(x > t) / rate_treated(x > t) <= target. Falls back to the
% strictest observed value when no threshold reaches the target.
    tr = ~is_ctrl;
    cand = sort(unique(x(tr & sig)), 'descend');
    thr = max(x) + 1;
    n_c = sum(is_ctrl); n_t = sum(tr);
    for t = cand(:)'
        above = x >= t & sig;
        fdr = (sum(above & is_ctrl) / n_c) / max(sum(above & tr) / n_t, eps);
        if fdr <= target
            thr = t - eps(t);   % calls are x > thr
        end
    end
    % note: keeps the least strict t that still meets the target (scan goes from strict to lax;
    % later qualifying t overwrite earlier ones)
end
