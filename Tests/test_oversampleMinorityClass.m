classdef test_oversampleMinorityClass < matlab.unittest.TestCase
% TEST_OVERSAMPLEMINORITYCLASS  Unit tests for Functions/oversampleMinorityClass.m
%
% Run: runtests('test_oversampleMinorityClass')

    methods (Test)

        function testImbalancedTwoClassOversampled(tc)
            % 4 class-1 + 2 class-2 training; 2 test rows
            n_total   = 8;
            n_train   = 6;
            feat      = array2table(rand(n_total, 3), 'VariableNames', {'f1','f2','f3'});
            grp       = 1:n_total;
            y_train   = [1 1 1 1 2 2];
            train_idx = logical([ones(1,n_train), zeros(1,n_total-n_train)]);

            [res_table, y_res, ~, tr_idx, te_idx] = ...
                oversampleMinorityClass(feat, grp, y_train, train_idx);

            % Both classes must reach the majority count (4)
            tc.verifyEqual(sum(y_res == 1), 4);
            tc.verifyEqual(sum(y_res == 2), 4);
            % Total = 8 train + 2 test
            tc.verifyEqual(height(res_table), 10);
            % Indices must be mutually exclusive and jointly exhaustive
            tc.verifyEqual(tr_idx & te_idx, false(1, height(res_table)));
            tc.verifyEqual(tr_idx | te_idx, true(1,  height(res_table)));
        end

        function testEqualClassesNoRowsAdded(tc)
            % 3+3 training already balanced — no oversampling expected
            n_each = 3; n_test = 2;
            n_total = 2*n_each + n_test;
            feat      = array2table(rand(n_total, 2), 'VariableNames', {'f1','f2'});
            grp       = 1:n_total;
            y_train   = [1 1 1 2 2 2];
            train_idx = logical([ones(1,2*n_each), zeros(1,n_test)]);

            [~, y_res, ~, tr_idx, ~] = ...
                oversampleMinorityClass(feat, grp, y_train, train_idx);

            tc.verifyEqual(sum(y_res == 1), n_each);
            tc.verifyEqual(sum(y_res == 2), n_each);
            tc.verifyEqual(sum(tr_idx), 2*n_each);
        end

        function testThreeClassImbalancedAllReachedMax(tc)
            % 6 / 3 / 2 → all should reach 6 after oversampling
            y_train   = [1 1 1 1 1 1  2 2 2  3 3];
            n_train   = length(y_train);
            n_test    = 2;
            n_total   = n_train + n_test;
            feat      = array2table(rand(n_total, 2), 'VariableNames', {'f1','f2'});
            grp       = 1:n_total;
            train_idx = logical([ones(1,n_train), zeros(1,n_test)]);

            [~, y_res, ~, ~, ~] = ...
                oversampleMinorityClass(feat, grp, y_train, train_idx);

            tc.verifyEqual(sum(y_res == 1), 6);
            tc.verifyEqual(sum(y_res == 2), 6);
            tc.verifyEqual(sum(y_res == 3), 6);
        end

        function testTestRowsUnchanged(tc)
            % Sentinel values in test rows must survive unchanged
            train_data = rand(4, 2);
            test_data  = [100, 200; 300, 400];
            feat       = array2table([train_data; test_data], 'VariableNames', {'f1','f2'});
            grp        = 1:6;
            y_train    = [1 1 2 2];
            train_idx  = logical([1 1 1 1 0 0]);

            [res_table, ~, ~, ~, te_idx] = ...
                oversampleMinorityClass(feat, grp, y_train, train_idx);

            test_rows = res_table(te_idx, :);
            tc.verifyEqual(test_rows.f1, [100; 300]);
            tc.verifyEqual(test_rows.f2, [200; 400]);
        end

        function testIndicesPartitionAllRows(tc)
            % train_idx and test_idx together must cover every row exactly once
            rng(11);
            feat      = array2table(rand(10, 3), 'VariableNames', {'f1','f2','f3'});
            grp       = 1:10;
            y_train   = [1 1 1 1 2 2 2];
            train_idx = logical([ones(1,7), zeros(1,3)]);

            [res_table, ~, ~, tr_idx, te_idx] = ...
                oversampleMinorityClass(feat, grp, y_train, train_idx);

            n = height(res_table);
            tc.verifyEqual(length(tr_idx), n);
            tc.verifyEqual(length(te_idx), n);
            tc.verifyEqual(tr_idx & te_idx, false(1, n));
            tc.verifyEqual(tr_idx | te_idx, true(1, n));
        end

        function testOutputTableHasOriginalColumns(tc)
            feat      = array2table(rand(6, 2), 'VariableNames', {'alpha','beta'});
            grp       = 1:6;
            y_train   = [1 1 2];
            train_idx = logical([1 1 1 0 0 0]);

            [res_table, ~, ~, ~, ~] = ...
                oversampleMinorityClass(feat, grp, y_train, train_idx);

            tc.verifyEqual(res_table.Properties.VariableNames, {'alpha','beta'});
        end

    end
end
