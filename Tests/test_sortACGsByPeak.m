classdef test_sortACGsByPeak < matlab.unittest.TestCase
% TEST_SORTACGSBYPEAK  Unit tests for Functions/sortACGsByPeak.m
%
% Run: runtests('test_sortACGsByPeak')

    methods (Test)

        function testOutputSizePreserved(tc)
            acgs = rand(5, 20);
            [sorted, idx] = sortACGsByPeak(acgs, "max");
            tc.verifySize(sorted, [5 20]);
            tc.verifySize(idx,    [5 1]);
        end

        function testSortedIdxIsValidPermutation(tc)
            rng(3);
            acgs = rand(8, 30);
            [sorted, idx] = sortACGsByPeak(acgs, "max");
            tc.verifyEqual(sort(idx), (1:8)', ...
                'sorted_idx must be a permutation of 1:N');
            tc.verifyEqual(sorted, acgs(idx, :));
        end

        function testDefaultEqualsMaxSort(tc)
            rng(5);
            acgs = rand(6, 40);
            [sorted_default, idx_default] = sortACGsByPeak(acgs);
            [sorted_max,     idx_max]     = sortACGsByPeak(acgs, "max");
            tc.verifyEqual(sorted_default, sorted_max);
            tc.verifyEqual(idx_default,    idx_max);
        end

        function testMeanSortOutputSize(tc)
            acgs = rand(7, 24);
            [sorted, idx] = sortACGsByPeak(acgs, "mean");
            tc.verifySize(sorted, [7 24]);
            tc.verifySize(idx,    [7 1]);
        end

        function testMeanSortIsValidPermutation(tc)
            rng(7);
            acgs = rand(5, 20);
            [~, idx] = sortACGsByPeak(acgs, "mean");
            tc.verifyEqual(sort(idx), (1:5)');
        end

        function testKnownOrderMaxSort(tc)
            % Row 1: peak at bin 2 (within first half of 10 bins → sort_key = 2)
            % Row 2: peak at bin 8 (outside first half   → sort_key = 1, argmax of zeros)
            % sort('descend') on [2;1] → row 1 comes first
            acgs = zeros(2, 10);
            acgs(1, 2) = 1;
            acgs(2, 8) = 1;
            [~, idx] = sortACGsByPeak(acgs, "max");
            tc.verifyEqual(idx(1), 1, ...
                'Row with peak in first half (sort_key=2) should rank first under descend');
        end

        function testSingleRowIdentity(tc)
            acgs = rand(1, 20);
            [sorted, idx] = sortACGsByPeak(acgs, "max");
            tc.verifyEqual(sorted, acgs);
            tc.verifyEqual(idx, 1);
        end

        function testUnknownTypeThrowsError(tc)
            acgs = rand(3, 10);
            tc.verifyError(@() sortACGsByPeak(acgs, "invalid_type"), '');
        end

        function testMaxSortDescendingByPeakBin(tc)
            % Construct ACGs where peak locations are 5, 3, 1 (in first half of 10)
            % sort('descend') should order rows as peak@5, peak@3, peak@1
            acgs = zeros(3, 10);
            acgs(1, 5) = 1;   % sort_key = 5
            acgs(2, 3) = 1;   % sort_key = 3
            acgs(3, 1) = 1;   % sort_key = 1
            [~, idx] = sortACGsByPeak(acgs, "max");
            tc.verifyEqual(idx, [1; 2; 3]);
        end

    end
end
