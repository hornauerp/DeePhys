classdef test_unitGraph < matlab.unittest.TestCase
% TEST_UNITGRAPH  Deterministic kNN graph and the community enrichment test.
%
% Run: runtests('test_unitGraph')

    methods (Test)

        function testGraphMatchesReferenceAndIsDeterministic(tc)
            rng(1);
            X = [randn(40, 3); randn(40, 3) + 10];
            k = 6;
            W1 = CellTypeClassifier.knnGraph(X, k);
            W2 = CellTypeClassifier.knnGraph(X(:, :), k);
            tc.verifyEqual(W1, W2);                               % deterministic
            tc.verifyEqual(full(W1), full(W1)', 'AbsTol', 1e-12); % symmetric
            tc.verifyEqual(full(diag(W1)), zeros(80, 1));
            tc.verifyLessThanOrEqual(full(max(W1(:))), 1 + 1e-12);
            tc.verifyEqual(full(W1), test_unitGraph.reference(X, k), 'AbsTol', 1e-6);
            % the two clusters are not connected at this k
            tc.verifyEqual(nnz(W1(1:40, 41:80)), 0);
        end

        function testDuplicateRowsHandled(tc)
            X = [zeros(3, 2); 1 0; 0 1; 1 1];
            W = CellTypeClassifier.knnGraph(X, 2);
            tc.verifyTrue(all(isfinite(nonzeros(W))));
            tc.verifyEqual(full(diag(W)), zeros(6, 1));
        end

        function testEnrichmentFindsResponderCommunity(tc)
            % three communities of 100 eligible units; responders concentrated in community 2
            M = repelem((1:3)', 100);
            resp = false(300, 1);
            resp(101:140) = true;            % 40 % in community 2
            resp([1:3, 201:202]) = true;     % background
            elig = true(300, 1);
            [T, inh, exc] = CellTypeClassifier.communityEnrichment(M, elig, resp, 0.05);
            tc.verifyEqual(inh, 2);
            tc.verifyTrue(all(ismember(exc, [1 3])));
            tc.verifyEqual(T.NResponsive', [3 40 2]);
            tc.verifyEqual(T.Expected(1), 300 * 45 / 300 / 3, 'AbsTol', 1e-12);
        end

        function testIneligibleUnitsIgnored(tc)
            % responders in ineligible units (other culture condition) do not count
            M = repelem((1:2)', 50);
            resp = false(100, 1); resp(51:70) = true;
            elig = true(100, 1); elig(51:100) = false;
            [T, inh] = CellTypeClassifier.communityEnrichment(M, elig, resp, 0.05);
            tc.verifyEqual(T.NResponsive', [0 0]);
            tc.verifyEmpty(inh);
        end

        function testSeedsBalancedWhenExcitatoryPoolIsSmall(tc)
            % pure-culture labels: 60 labelled inhibitory units, only 20 labelled excitatory
            % ones. The propagation compares per-seed absorption, so the seed classes must be
            % equal in size; the surplus inhibitory seeds are dropped.
            rng(5);
            X = [randn(100, 2) + [4 0]; randn(100, 2) - [4 0]];      % I cluster, E cluster
            G = CellTypeClassifier.knnGraph(X, 10);
            resp = false(200, 1); resp(1:60) = true;                  % labelled I
            elig = resp; elig(101:120) = true;                        % plus 20 labelled E
            p_comm = struct('LouvainResolution', 1, 'LouvainRestarts', 10, 'SeedFDR', 0.05, 'Algorithm', "leiden");
            [inh, exc] = CellTypeClassifier.enrichmentSeeds(G, resp, elig, p_comm, 1);
            tc.verifyEqual(numel(inh), numel(exc));
            tc.verifyTrue(all(inh <= 60));
            tc.verifyTrue(all(exc >= 101 & exc <= 120));
        end

        function testNoStructureGivesNoEnrichment(tc)
            rng(3);
            M = repelem((1:5)', 60);
            resp = rand(300, 1) < 0.1;
            [~, inh] = CellTypeClassifier.communityEnrichment(M, true(300, 1), resp, 0.05);
            tc.verifyEmpty(inh);
        end

    end

    methods (Test)

        function testWNNFollowsInformativeModality(tc)
            % Units of three types: the "ACG" block separates them, the "waveform" block is
            % pure high-variance noise. Fixed equal weighting lets the noise dominate; WNN
            % gives the ACG block most of the weight and recovers the types.
            rng(31);
            n = 300; type = repelem((1:3)', 100);
            centers = 3 * randn(3, 30);
            acg = centers(type, :) + randn(n, 30);
            wf  = 15 * randn(n, 30);
            X = [acg, wf]; is_acg = [true(1, 30), false(1, 30)];
            [W, info] = CellTypeClassifier.wnnGraph(X, is_acg, 15);
            tc.verifyEqual(W, W', 'AbsTol', 1e-12);
            tc.verifyEqual(full(diag(W)), zeros(n, 1));
            tc.verifyTrue(all(nonzeros(W) > 0 & nonzeros(W) <= 1));
            tc.verifyGreaterThan(median(info.weights(:, 1)), 0.6);
            purity = @(G) full(sum(sum(G .* (type == type')))) / full(sum(G(:)));
            G_cat = CellTypeClassifier.knnGraph(X, 15);   % concatenated features
            tc.verifyGreaterThan(purity(W), 0.9);
            tc.verifyGreaterThan(purity(W), purity(G_cat) + 0.2);
        end

        function testWNNPerUnitWeights(tc)
            % First half: types visible in modality A only; second half: in modality B only.
            rng(32);
            n = 400; type = repmat(repelem((1:2)', 100), 2, 1);
            first = (1:n)' <= 200;
            ca = 4 * randn(2, 20); cb = 4 * randn(2, 20);
            A = randn(n, 20) * 2; B = randn(n, 20) * 2;
            A(first, :)  = ca(type(first), :) + randn(200, 20);
            B(~first, :) = cb(type(~first), :) + randn(200, 20);
            [~, info] = CellTypeClassifier.wnnGraph([A, B], [true(1, 20), false(1, 20)], 15);
            tc.verifyGreaterThan(median(info.weights(first, 1)), 0.5);
            tc.verifyLessThan(median(info.weights(~first, 1)), 0.5);
        end

    end

    methods (Static)

        function W = reference(X, k)
            % Straightforward UMAP fuzzy simplicial set with exact neighbours (loops)
            N = size(X, 1); D = squareform(pdist(X));
            A = zeros(N);
            for i = 1:N
                d = D(i, :); d(i) = Inf;
                [ds, o] = sort(d); ds = ds(1:k); o = o(1:k);
                rho = min(ds(ds > 0));
                f = @(s) sum(exp(-max(ds - rho, 0) / s)) - log2(k);
                s = fzero(f, [1e-6, 1e3]);
                A(i, o) = exp(-max(ds - rho, 0) / s);
            end
            W = A + A' - A .* A';
        end

    end
end
