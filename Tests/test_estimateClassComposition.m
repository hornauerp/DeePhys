classdef test_estimateClassComposition < matlab.unittest.TestCase
% TEST_ESTIMATECLASSCOMPOSITION  Group fractions under prior shift, with unit identity kept.
%
% Run: runtests('test_estimateClassComposition')

    methods (Test)

        function testRecoversGroupFractionsFromPureTraining(tc)
            % two Gaussian classes in 5 D; training only on pure groups (as pure cultures),
            % applied to groups of 10 / 50 / 80 % positives
            rng(31);
            truth = [0 1 0.1 0.5 0.8];
            Z = []; grp = strings(0, 1); lab = [];
            for g = 1:numel(truth)
                n = 400; yy = rand(n, 1) < truth(g);
                Z = [Z; randn(n, 5) + 1.2 * (2 * yy - 1) .* [1 0.5 0 0 0]]; %#ok<AGROW>
                grp = [grp; repmat("g" + g, n, 1)]; lab = [lab; yy]; %#ok<AGROW>
            end
            train = find(ismember(grp, ["g1", "g2"]));
            [T, post, p] = estimateClassComposition(Z, train, 1 + lab(train), grp, NBoot=100);
            est = T.Fraction(ismember(T.Group, ["g3", "g4", "g5"]));
            tc.verifyEqual(est', truth(3:5), 'AbsTol', 0.07);
            tc.verifyTrue(all(T.CILow <= T.Fraction & T.Fraction <= T.CIHigh));
            m = ismember(grp, ["g3", "g4", "g5"]);
            auc = @(s, y) (sum(tiedrank(s) .* y) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(~y));
            tc.verifyGreaterThan(auc(post(m), lab(m)), 0.85);      % unit identity is kept
            % the prior adjustment beats counting at the training prior in the extreme group
            g3 = grp == "g3";
            tc.verifyLessThan(abs(T.Fraction(T.Group == "g3") - 0.1), abs(mean(p(g3) > 0.5) - 0.1));
        end

        function testModelScoresNewRows(tc)
            % the returned model gives the probabilities of the non-training rows
            rng(32);
            y = rand(300, 1) < 0.5; Z = randn(300, 3) + 1.5 * y;
            grp = repmat("g", 300, 1);
            [~, ~, p, mdl] = estimateClassComposition(Z, 1:150, 1 + y(1:150), grp, NBoot=5);
            [~, s] = predict(mdl, Z(151:end, :));
            tc.verifyEqual(s(:, 2), p(151:end), 'AbsTol', 1e-12);
        end

        function testNeedsBothClasses(tc)
            tc.verifyError(@() estimateClassComposition(randn(10, 2), 1:5, ones(5, 1), ...
                repmat("a", 10, 1)), 'estimateClassComposition:oneClass');
        end

    end
end
