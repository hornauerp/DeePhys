function [T, inh_comm, exc_comm] = communityEnrichmentTest(M, elig, resp, fdr)
% COMMUNITYENRICHMENTTEST  Hypergeometric responder enrichment / depletion per community.
%
%   [T, inh_comm, exc_comm] = communityEnrichmentTest(M, elig, resp, fdr)
%
% INPUTS:
%   M    - (N x 1) community index per unit (1..K)
%   elig - (N x 1) logical, unit is eligible (in a culture where responders occur)
%   resp - (N x 1) logical, unit is responsive (only counted when eligible)
%   fdr  - Benjamini-Hochberg false-discovery rate across communities
% OUTPUTS:
%   T        - table: Community, NEligible, NResponsive, Expected, pEnriched, pDepleted,
%              qEnriched, qDepleted
%   inh_comm - communities significantly enriched in responders (qEnriched < fdr)
%   exc_comm - communities significantly depleted (qDepleted < fdr)
    M = M(:); elig = logical(elig(:)); resp = logical(resp(:)) & elig;
    K = max(M);
    N_e = sum(elig); R_e = sum(resp);
    n_c = accumarray(M(elig), 1, [K 1]);
    r_c = accumarray(M(resp), 1, [K 1]);
    tested = n_c > 0;
    p_enr = ones(K, 1); p_dep = ones(K, 1);
    p_enr(tested) = hygecdf(r_c(tested) - 1, N_e, R_e, n_c(tested), 'upper');
    p_dep(tested) = hygecdf(r_c(tested), N_e, R_e, n_c(tested));
    q_enr = ones(K, 1); q_dep = ones(K, 1);
    q_enr(tested) = bhAdjust(p_enr(tested));
    q_dep(tested) = bhAdjust(p_dep(tested));
    T = table((1:K)', n_c, r_c, n_c * R_e / max(N_e, 1), p_enr, p_dep, q_enr, q_dep, ...
        'VariableNames', {'Community', 'NEligible', 'NResponsive', 'Expected', ...
                          'pEnriched', 'pDepleted', 'qEnriched', 'qDepleted'});
    inh_comm = find(tested & q_enr < fdr)';
    exc_comm = find(tested & q_dep < fdr)';
end

function q = bhAdjust(p)
% Benjamini-Hochberg adjusted p-values (column in, column out)
    m = numel(p);
    [ps, o] = sort(p(:));
    qs = ps .* m ./ (1:m)';
    qs = flipud(cummin(flipud(qs)));
    q = zeros(m, 1); q(o) = min(qs, 1);
end
