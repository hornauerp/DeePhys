classdef test_RegressorPermutation < matlab.unittest.TestCase
% TEST_REGRESSORPERMUTATION  Regression permutation test respects culture groups.
%
% Run: runtests('test_RegressorPermutation')

    methods (Test)

        function testCultureConfoundNotSignificant(tc)
            % Target constant per culture (e.g. age), features carry only culture identity:
            % shuffling single units would call this significant, whole cultures must not.
            rng(21);
            [X, Y, cult] = test_RegressorPermutation.makeCultures(8, 30, 0);
            opts = struct('Algorithm', 'knn', 'KFold', 4, 'Seed', 1, 'CVGroups', cult);
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            p = Regressor.permutationTest(X, Y, opts, 19);
            tc.verifyGreaterThan(p, 0.05);
        end

        function testPValueNeverZero(tc)
            rng(22);
            x = randn(80, 1);
            X = table(x); Y = 3 * x + 0.1 * randn(80, 1);
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            p = Regressor.permutationTest(X, Y, struct('Algorithm', 'knn', 'Seed', 1), 9);
            tc.verifyEqual(p, 0.1, 'AbsTol', 1e-12);
        end

        function testLabelPermuterSchemes(tc)
            g = repelem((1:4)', 3);
            [f, s] = MLUtils.labelPermuter(repelem([1; 2; 3; 4], 3), g);
            tc.verifyEqual(s, "group");
            y = f();
            tc.verifyEqual(sort(unique(y)), (1:4)');
            tc.verifyTrue(all(splitapply(@(v) all(v == v(1)), y, g)));   % still constant per group
            [f, s] = MLUtils.labelPermuter((1:12)', g);
            tc.verifyEqual(s, "within-group");
            y = f();
            tc.verifyEqual(splitapply(@(v) {sort(v)}, y, g), splitapply(@(v) {sort(v)}, (1:12)', g));
            [~, s] = MLUtils.labelPermuter((1:12)', []);
            tc.verifyEqual(s, "observation");
        end

    end

    methods (Static)

        function [X, Y, cult] = makeCultures(n_cult, n_units, effect)
            cult = repelem("c" + (1:n_cult)', n_units);
            age  = repelem((1:n_cult)', n_units);           % culture-level target
            offs = repelem(2 * randn(n_cult, 3), n_units, 1);
            F = randn(n_cult * n_units, 3) * 0.3 + offs;
            F(:, 1) = F(:, 1) + effect * age;
            X = array2table(F, 'VariableNames', ["f1", "f2", "f3"]);
            Y = age;
        end

    end
end
