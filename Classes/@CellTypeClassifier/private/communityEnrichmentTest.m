function [T, inh_comm, exc_comm] = communityEnrichmentTest(M, elig, resp, fdr, strata)
% COMMUNITYENRICHMENTTEST  Responder enrichment / depletion per community, stratified by culture.
%
%   [T, inh_comm, exc_comm] = communityEnrichmentTest(M, elig, resp, fdr, strata)
%
% Per community, the responder count among its eligible units is tested against the exact
% conditional null of a stratified (Cochran-Mantel-Haenszel-type) test: within each stratum
% (culture) the count is hypergeometric given the stratum's responders and the community's
% units in it, and the strata are summed by convolution. Pooling cultures instead would call
% a community "enriched" merely because it holds many units of a culture with a high responder
% rate (expression level, drift). With a single stratum this is the plain hypergeometric test.
%
% INPUTS:
%   M      - (N x 1) community index per unit (1..K)
%   elig   - (N x 1) logical, unit is eligible (in a culture where responders occur)
%   resp   - (N x 1) logical, unit is responsive (only counted when eligible)
%   fdr    - Benjamini-Hochberg false-discovery rate across communities
%   strata - (N x 1) stratum (culture) per unit; [] = one stratum
% OUTPUTS:
%   T        - table: Community, NEligible, NResponsive, Expected, pEnriched, pDepleted,
%              qEnriched, qDepleted
%   inh_comm - communities significantly enriched in responders (qEnriched < fdr)
%   exc_comm - communities significantly depleted (qDepleted < fdr)
    if nargin < 5 || isempty(strata), strata = ones(numel(M), 1); end
    M = M(:); elig = logical(elig(:)); resp = logical(resp(:)) & elig;
    [~, ~, s] = unique(string(strata(:)));
    K = max(M);
    S = max([s(elig); 0]);
    % per stratum: eligible units and responders; per community x stratum: units and responders
    N_s = accumarray(s(elig), 1, [S 1]);
    R_s = accumarray(s(resp), 1, [S 1]);
    n_cs = accumarray([M(elig) s(elig)], 1, [K S]);
    r_c  = accumarray(M(resp), 1, [K 1]);
    n_c  = sum(n_cs, 2);
    expected = n_cs * (R_s ./ max(N_s, 1));
    tested = n_c > 0;
    p_enr = ones(K, 1); p_dep = ones(K, 1);
    for c = find(tested)'
        [p_enr(c), p_dep(c)] = stratifiedHypergeomTails(r_c(c), N_s, R_s, n_cs(c, :)');
    end
    q_enr = ones(K, 1); q_dep = ones(K, 1);
    q_enr(tested) = bhAdjust(p_enr(tested));
    q_dep(tested) = bhAdjust(p_dep(tested));
    T = table((1:K)', n_c, r_c, expected, p_enr, p_dep, q_enr, q_dep, ...
        'VariableNames', {'Community', 'NEligible', 'NResponsive', 'Expected', ...
                          'pEnriched', 'pDepleted', 'qEnriched', 'qDepleted'});
    inh_comm = find(tested & q_enr < fdr)';
    exc_comm = find(tested & q_dep < fdr)';
end

function [p_upper, p_lower] = stratifiedHypergeomTails(a, N, R, n)
% P(A >= a) and P(A <= a) for A = sum of independent hypergeometric counts, one per stratum
% (population N(k), R(k) successes, n(k) draws).
    pmf = 1; lo = 0;
    for k = find(n(:)' > 0)
        x = max(0, n(k) + R(k) - N(k)):min(n(k), R(k));
        pmf = conv(pmf, hygepdf(x, N(k), R(k), n(k)));
        lo = lo + x(1);
    end
    idx = a - lo + 1;                       % position of a in the support
    if idx < 1
        p_upper = 1; p_lower = 0;
    elseif idx > numel(pmf)
        p_upper = 0; p_lower = 1;
    else
        p_upper = min(1, sum(pmf(idx:end)));
        p_lower = min(1, sum(pmf(1:idx)));
    end
end

function q = bhAdjust(p)
% Benjamini-Hochberg adjusted p-values (column in, column out)
    m = numel(p);
    [ps, o] = sort(p(:));
    qs = ps .* m ./ (1:m)';
    qs = flipud(cummin(flipud(qs)));
    q = zeros(m, 1); q(o) = min(qs, 1);
end
