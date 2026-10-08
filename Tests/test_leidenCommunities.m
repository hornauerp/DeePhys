classdef test_leidenCommunities < matlab.unittest.TestCase
% TEST_LEIDENCOMMUNITIES  Leiden modularity optimisation (leidenCommunities).
%
% Run: runtests('test_leidenCommunities')

    methods (Test)

        function testKarateClubReachesKnownOptimum(tc)
            A = test_leidenCommunities.karate();
            Q = zeros(1, 10);
            for r = 1:10
                rng(r);
                [M, Q(r)] = leidenCommunities(A);
                tc.verifyEqual(Q(r), test_leidenCommunities.modularity(A, M, 1), 'AbsTol', 1e-12);
                tc.verifyTrue(test_leidenCommunities.allConnected(A, M));
            end
            tc.verifyGreaterThan(max(Q), 0.4197);   % optimum 0.41979 (4 communities)
        end

        function testRecoversPlantedBlocks(tc)
            rng(4);
            n = 200; truth = repelem((1:4)', 50);
            P = 0.02 + 0.28 * (truth == truth');
            A = triu(rand(n) < P, 1) .* rand(n);  % weighted
            A = sparse(A + A');
            [M, Q] = leidenCommunities(A);
            tc.verifyEqual(max(M), 4);
            tc.verifyEqual(test_leidenCommunities.ari(M, truth), 1, 'AbsTol', 1e-12);
            tc.verifyTrue(test_leidenCommunities.allConnected(A, M));
            tc.verifyEqual(Q, test_leidenCommunities.modularity(A, M, 1), 'AbsTol', 1e-12);
        end

        function testResolutionControlsGranularity(tc)
            rng(5);
            A = test_leidenCommunities.karate();
            M_lo = leidenCommunities(A, 0.5);
            M_hi = leidenCommunities(A, 2);
            tc.verifyLessThan(max(M_lo), max(M_hi));
            [M0, Q0] = leidenCommunities(A, 1, Initial=M_hi);   % refine from a given partition
            tc.verifyGreaterThanOrEqual(Q0, test_leidenCommunities.modularity(A, M_hi, 1) - 1e-12);
            tc.verifyEqual(numel(M0), 34);
        end

        function testEmptyGraph(tc)
            [M, Q] = leidenCommunities(sparse(5, 5));
            tc.verifyEqual(M, (1:5)');
            tc.verifyEqual(Q, 0);
        end

    end

    methods (Static)

        function A = karate()
            % Zachary's karate club (34 members, 78 ties)
            E = [1 2;1 3;1 4;1 5;1 6;1 7;1 8;1 9;1 11;1 12;1 13;1 14;1 18;1 20;1 22;1 32; ...
                 2 3;2 4;2 8;2 14;2 18;2 20;2 22;2 31;3 4;3 8;3 9;3 10;3 14;3 28;3 29;3 33; ...
                 4 8;4 13;4 14;5 7;5 11;6 7;6 11;6 17;7 17;9 31;9 33;9 34;10 34;14 34; ...
                 15 33;15 34;16 33;16 34;19 33;19 34;20 34;21 33;21 34;23 33;23 34; ...
                 24 26;24 28;24 30;24 33;24 34;25 26;25 28;25 32;26 32;27 30;27 34;28 34; ...
                 29 32;29 34;30 33;30 34;31 33;31 34;32 33;32 34;33 34];
            A = sparse(E(:, 1), E(:, 2), 1, 34, 34);
            A = A + A';
        end

        function Q = modularity(A, M, gamma)
            k = full(sum(A, 2)); m2 = sum(k);
            same = M(:) == M(:)';
            Q = full(sum(sum((A - gamma * (k * k') / m2) .* same))) / m2;
        end

        function tf = allConnected(A, M)
            tf = true;
            for c = 1:max(M)
                idx = find(M == c);
                if numel(idx) > 1
                    bins = conncomp(graph(A(idx, idx) ~= 0));
                    tf = tf && max(bins) == 1;
                end
            end
        end

        function a = ari(p, q)
            [~, ~, p] = unique(p(:)); [~, ~, q] = unique(q(:));
            N = accumarray([p q], 1); c2 = @(x) x .* (x - 1) / 2;
            s = sum(c2(N(:))); sa = sum(c2(sum(N, 2))); sb = sum(c2(sum(N, 1)));
            e = sa * sb / c2(numel(p));
            a = (s - e) / ((sa + sb) / 2 - e);
        end

    end
end
