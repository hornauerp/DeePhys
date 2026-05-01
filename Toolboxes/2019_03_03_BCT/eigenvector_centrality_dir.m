function v = eigenvector_centrality_dir(CIJ)
%EIGENVECTOR_CENTRALITY_DIR      Spectral measure of centrality for directed graphs
%
%   v = eigenvector_centrality_dir(CIJ)
%
%   Eigenvector centrality is a self-referential measure of centrality:
%   nodes have high eigenvector centrality if they connect to other nodes
%   that have high eigenvector centrality. The eigenvector centrality of
%   node i is equivalent to the ith element in the eigenvector
%   corresponding to the largest eigenvalue of the adjacency matrix.
%
%   For directed (asymmetric) networks, eigenvalues may be complex.
%   The eigenvector associated with the eigenvalue of largest magnitude
%   is selected.
%
%   Inputs:     CIJ,        binary/weighted directed adjacency matrix
%                           (possibly asymmetric).
%
%   Outputs:      v,        eigenvector centrality (n x 1, real non-negative).
%
%   Note: For undirected graphs, use eigenvector_centrality_und instead.
%
%   Reference: Newman, MEJ (2002). The mathematics of networks.
%              Bonacich, P (1972). J Math Soc, 2(1), 113-120.

n = length(CIJ);
if n < 1000
    [V, D] = eig(CIJ);
else
    [V, D] = eigs(sparse(CIJ));
end
[~, idx] = max(abs(diag(D)));
ec = abs(V(:, idx));
v = reshape(ec, length(ec), 1);
