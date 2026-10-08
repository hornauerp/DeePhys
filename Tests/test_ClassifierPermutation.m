classdef test_ClassifierPermutation < matlab.unittest.TestCase
% TEST_CLASSIFIERPERMUTATION  Permutation test respects the group structure of labels.
%
% Run: runtests('test_ClassifierPermutation')

    methods (Test)

        function testCultureConfoundNotSignificant(tc)
            % Labels are constant per culture and features carry only culture identity
            % (no label signal). Group-level permutation must not call this significant.
            rng(11);
            [X, Y, cult] = test_ClassifierPermutation.makeCultures(8, 40, 0);
            opts = struct('Algorithm', 'knn', 'KFold', 4, 'Seed', 1, 'CVGroups', cult);
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            p = Classifier.permutationTest(X, Y, opts, 19);
            tc.verifyGreaterThan(p, 0.05);
        end

        function testRealSignalIsSignificant(tc)
            % 12 cultures: a group permutation reproduces the true split (or its mirror)
            % with probability 2/924, so the null rarely contains the observed labelling.
            rng(12);
            [X, Y, cult] = test_ClassifierPermutation.makeCultures(12, 30, 10);
            opts = struct('Algorithm', 'knn', 'KFold', 4, 'Seed', 1, 'CVGroups', cult);
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            [p, ~, obs] = Classifier.permutationTest(X, Y, opts, 19);
            tc.verifyLessThanOrEqual(p, 0.1);
            tc.verifyGreaterThan(obs, 0.8);
        end

        function testPValueNeverZero(tc)
            rng(13);
            X = table([randn(30, 1); randn(30, 1) + 5]);
            Y = [zeros(30, 1); ones(30, 1)];
            w = warning('off', 'all');
            c = onCleanup(@() warning(w));
            p = Classifier.permutationTest(X, Y, struct('Algorithm', 'knn', 'Seed', 1), 9);
            tc.verifyEqual(p, 0.1, 'AbsTol', 1e-12);
        end

    end

    methods (Static)

        function [X, Y, cult] = makeCultures(n_cult, n_units, effect)
            % Half the cultures are "I", half "E". Each culture has a strong random offset
            % (batch effect); "effect" adds a true label difference on f1.
            cult = repelem("c" + (1:n_cult)', n_units);
            Y = repelem(repmat(["E"; "I"], n_cult / 2, 1), n_units);
            offs = repelem(2 * randn(n_cult, 3), n_units, 1);
            F = randn(n_cult * n_units, 3) * 0.3 + offs;
            F(:, 1) = F(:, 1) + effect * (Y == "I");
            X = array2table(F, 'VariableNames', ["f1", "f2", "f3"]);
        end

    end
end
