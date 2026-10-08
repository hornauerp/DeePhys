classdef test_zscoreByGroup < matlab.unittest.TestCase
% TEST_ZSCOREBYGROUP  Per-group z-score, plain and with bridge reference rows.
%
% Run: runtests('test_zscoreByGroup')

    methods (Test)

        function testPlainMatchesNormalize(tc)
            rng(1);
            X = randn(60, 3) + [0 5 -2]; g = repelem(["a"; "b"; "c"], 20);
            Z = zscoreByGroup(X, g);
            for k = ["a", "b", "c"]
                tc.verifyEqual(Z(g == k, :), normalize(X(g == k, :)), 'AbsTol', 1e-12);
            end
        end

        function testBridgeKeepsCompositionDifference(tc)
            % Two plates with a technical offset of +3. Both carry bridge wells of equal
            % composition; plate 1 also has wells with more "high" units. Plain per-plate
            % z-scoring removes the composition difference, the bridge reference keeps it.
            rng(2);
            bridge = [randn(100, 1); randn(100, 1) + 2];             % 50:50 low/high
            rich   = [randn(30, 1);  randn(170, 1) + 2];             % mostly high
            poor   = [randn(170, 1); randn(30, 1) + 2];              % mostly low
            X = [bridge; rich; bridge + 3; poor + 3];
            plate  = [ones(400, 1); 2 * ones(400, 1)];
            is_ref = [true(200, 1); false(200, 1); true(200, 1); false(200, 1)];
            Zb = zscoreByGroup(X, plate, is_ref);
            Zp = zscoreByGroup(X, plate);
            gap_bridge = mean(Zb(201:400)) - mean(Zb(601:800));
            gap_plain  = mean(Zp(201:400)) - mean(Zp(601:800));
            tc.verifyGreaterThan(gap_bridge, 0.75);                   % composition difference kept (expected 1.4 / sqrt(2))
            tc.verifyLessThan(abs(mean(Zb(1:200)) - mean(Zb(401:600))), 1e-10);   % offset removed
            tc.verifyGreaterThan(gap_bridge, 1.5 * gap_plain);
        end

        function testFallbackWithoutReferenceRows(tc)
            X = [randn(30, 2); randn(30, 2) + 4];
            g = [ones(30, 1); 2 * ones(30, 1)];
            is_ref = [true(30, 1); false(30, 1)];                      % group 2 has no reference
            [Z, fb] = zscoreByGroup(X, g, is_ref);
            tc.verifyEqual(fb, 2);
            tc.verifyEqual(Z(31:60, :), normalize(X(31:60, :)), 'AbsTol', 1e-12);
        end

        function testStoredStatsStandardiseNewRows(tc)
            % rows standardised with the stored statistics equal the fitted ones; unknown
            % groups become NaN
            rng(3);
            X = randn(60, 3) + [0 5 -2]; g = repelem(["a"; "b"; "c"], 20);
            is_ref = rand(60, 1) < 0.6;
            [Z, ~, S] = zscoreByGroup(X, g, is_ref, 5);
            Z2 = zscoreByGroup(X(41:60, :), g(41:60), Stats=S);
            tc.verifyEqual(Z2, Z(41:60, :), 'AbsTol', 1e-12);
            tc.verifyTrue(all(isnan(zscoreByGroup(X(1:2, :), ["z"; "z"], Stats=S)), 'all'));
        end

        function testConstantColumnBecomesZero(tc)
            X = [ones(25, 1), (1:25)'];
            Z = zscoreByGroup(X, ones(25, 1));
            tc.verifyEqual(Z(:, 1), zeros(25, 1));
            % ConstantAs NaN reproduces MATLAB's normalize, group by group
            X2 = [X; X + [0 3]]; g = [ones(25, 1); 2 * ones(25, 1)];
            N = zscoreByGroup(X2, g, true(50, 1), 20, ConstantAs=NaN);
            tc.verifyEqual(N, [normalize(X2(1:25, :)); normalize(X2(26:50, :))], 'AbsTol', 1e-12);
        end

    end
end
