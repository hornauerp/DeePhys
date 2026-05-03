classdef test_GroupedCV < matlab.unittest.TestCase
% TEST_GROUPEDCV  Unit tests for GroupedCV hierarchy-aware cross-validation.
%
% Run: runtests('test_GroupedCV')

    methods (Test)

        function testFoldMasksAreComplementary(tc)
            group_ids = repelem([1;2;3;4;5], 10);  % 5 groups of 10
            gcv = GroupedCV(group_ids, 5);
            for k = 1:gcv.NumFolds
                [tr, te] = gcv.fold(k);
                tc.verifyEqual(tr | te, true(50, 1));
                tc.verifyEqual(tr & te, false(50, 1));
            end
        end

        function testNoGroupLeak(tc)
            % Units from the same recording must not appear in both train and test
            group_ids = repelem([1;2;3;4;5], 10);
            gcv = GroupedCV(group_ids, 5);
            for k = 1:gcv.NumFolds
                [tr, te] = gcv.fold(k);
                groups_in_train = unique(group_ids(tr));
                groups_in_test  = unique(group_ids(te));
                tc.verifyEmpty(intersect(groups_in_train, groups_in_test));
            end
        end

        function testLeaveOneGroupOut(tc)
            group_ids = repelem((1:8)', 5);  % 8 groups of 5
            gcv = GroupedCV(group_ids, -1);
            tc.verifyEqual(gcv.NumFolds, 8);
            for k = 1:gcv.NumFolds
                [~, te] = gcv.fold(k);
                tc.verifyEqual(sum(te), 5);  % exactly one group per test fold
            end
        end

        function testReducesToNGroupsWhenKTooLarge(tc)
            group_ids = repelem((1:3)', 10);
            % Request 10 folds with only 3 groups — should warn and reduce
            tc.verifyWarning( ...
                @() GroupedCV(group_ids, 10), ...
                'GroupedCV:foldReduction');
            warnState = warning('off', 'GroupedCV:foldReduction');
            gcv = GroupedCV(group_ids, 10);
            warning(warnState);
            tc.verifyEqual(gcv.NumFolds, 3);
        end

        function testStratificationPreservesClassBalance(tc)
            % 4 groups: 2 class-A, 2 class-B → each fold should have both classes in test
            group_ids = repelem([1;2;3;4], 10);
            Y = [ones(20,1); 2*ones(20,1)];
            gcv = GroupedCV(group_ids, 2, Y);
            tc.verifyEqual(gcv.NumFolds, 2);
            for k = 1:gcv.NumFolds
                [~, te] = gcv.fold(k);
                y_te = Y(te);
                tc.verifyTrue(any(y_te == 1));
                tc.verifyTrue(any(y_te == 2));
            end
        end

        function testByGroupsFactoryWithStrings(tc)
            group_ids = repelem(["rec1","rec2","rec3","rec4"], 10)';
            Y = repelem([1;2;1;2], 10);
            gcv = GroupedCV.byGroups(group_ids, 4, Y);
            tc.verifyEqual(gcv.NumFolds, 4);
        end

        function testDescribeSummary(tc)
            group_ids = repelem((1:5)', 8);
            gcv = GroupedCV(group_ids, 5);
            s = gcv.describe();
            tc.verifyEqual(s.NumFolds, 5);
            tc.verifyEqual(s.NumGroups, 5);
            tc.verifyEqual(s.NumObjects, 40);
            tc.verifyEqual(s.MeanGroupSize, 8);
        end

        function testTwoGroupsWithKTooLargeReducesTo2(tc)
            group_ids = repelem([1;2], 10);
            warnState = warning('off', 'GroupedCV:foldReduction');
            gcv = GroupedCV(group_ids, 10);
            warning(warnState);
            tc.verifyEqual(gcv.NumFolds, 2);
            % Both folds should still work
            for k = 1:gcv.NumFolds
                [tr, te] = gcv.fold(k);
                tc.verifyTrue(sum(tr) > 0 && sum(te) > 0);
            end
        end

    end
end
