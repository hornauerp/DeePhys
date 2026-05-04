classdef test_CellTypeClassifier_helpers < matlab.unittest.TestCase
% TEST_CELLTYPECLASSIFIER_HELPERS  Unit tests for CellTypeClassifier static helpers.
%
% Tests the three static methods that have no external toolbox dependencies:
%   evaluateLabels   — classification scoring (pure computation)
%   computeACG       — autocorrelogram via CCG.m (requires Functions on path)
%   binnedSpikeMatrix — spike binning from UnitData objects
%
% Run: runtests('test_CellTypeClassifier_helpers')

    methods (Test)

        % ------------------------------------------------------------------
        % evaluateLabels
        % ------------------------------------------------------------------

        function testEvaluatePerfectClassification(tc)
            stats = CellTypeClassifier.evaluateLabels([1 1 2 2], [1 1 2 2]);
            tc.verifyEqual(stats.accuracy, 1.0);
            tc.verifyEqual(stats.confusionMatrix, [2 0; 0 2]);
            tc.verifyEqual(stats.perClassAccuracy, [1.0, 1.0]);
            tc.verifyEqual(stats.nValid, 4);
            tc.verifyEqual(stats.nExcluded, 0);
        end

        function testEvaluateAllWrong(tc)
            stats = CellTypeClassifier.evaluateLabels([2 2 1 1], [1 1 2 2]);
            tc.verifyEqual(stats.accuracy, 0.0);
            tc.verifyEqual(stats.confusionMatrix, [0 2; 2 0]);
        end

        function testEvaluateNaNExclusion(tc)
            stats = CellTypeClassifier.evaluateLabels([1 NaN 2], [1 1 2]);
            tc.verifyEqual(stats.nExcluded, 1);
            tc.verifyEqual(stats.nValid, 2);
            tc.verifyEqual(stats.accuracy, 1.0);
        end

        function testEvaluateMismatchedLengthErrors(tc)
            tc.verifyError( ...
                @() CellTypeClassifier.evaluateLabels([1 2], [1 2 1]), ...
                'MATLAB:assert:failed');
        end

        function testEvaluateAllNaNReturnsZero(tc)
            stats = CellTypeClassifier.evaluateLabels([NaN NaN], [1 2]);
            tc.verifyEqual(stats.nValid, 0);
            tc.verifyEqual(stats.accuracy, 0.0);
        end

        % ------------------------------------------------------------------
        % computeACG
        % ------------------------------------------------------------------

        function testACGOutputSize(tc)
            rng(1);
            spike_times = sort(rand(500, 1) * 60);
            bin_size = 0.0005;  lag = 0.1;
            acg = CellTypeClassifier.computeACG(spike_times, bin_size, lag);
            expected_len = round(lag / bin_size) * 2 + 1;
            tc.verifyEqual(numel(acg), expected_len);
        end

        function testACGCenterBinIsZero(tc)
            rng(2);
            spike_times = sort(rand(500, 1) * 60);
            bin_size = 0.0005;  lag = 0.1;
            acg = CellTypeClassifier.computeACG(spike_times, bin_size, lag);
            center = (numel(acg) + 1) / 2;
            tc.verifyEqual(acg(center), 0);
        end

        function testACGEmptySpikesReturnsZeros(tc)
            bin_size = 0.0005;  lag = 0.1;
            acg_empty = CellTypeClassifier.computeACG(zeros(0,1), bin_size, lag);
            acg_one   = CellTypeClassifier.computeACG([0.5], bin_size, lag);
            expected_len = round(lag / bin_size) * 2 + 1;
            tc.verifyEqual(numel(acg_empty), expected_len);
            tc.verifyEqual(sum(acg_empty), 0);
            tc.verifyEqual(numel(acg_one), expected_len);
            tc.verifyEqual(sum(acg_one), 0);
        end

        function testACGCustomParams(tc)
            rng(3);
            spike_times = sort(rand(200, 1) * 60);
            bin_size = 0.001;  lag = 0.05;
            acg = CellTypeClassifier.computeACG(spike_times, bin_size, lag);
            expected_len = round(lag / bin_size) * 2 + 1;
            tc.verifyEqual(numel(acg), expected_len);
        end

        % ------------------------------------------------------------------
        % binnedSpikeMatrix
        % ------------------------------------------------------------------

        function testBinnedMatrixOutputDimensions(tc)
            ud = test_CellTypeClassifier_helpers.makeUnitData([0.5, 1.5, 2.5]);
            ud(2) = test_CellTypeClassifier_helpers.makeUnitData([0.1, 3.5]);
            ud(3) = test_CellTypeClassifier_helpers.makeUnitData([]);
            binned = CellTypeClassifier.binnedSpikeMatrix(ud, [0, 5], 1.0);
            tc.verifySize(binned, [3, 5]);
        end

        function testBinnedMatrixSpikeCountsCorrect(tc)
            % One unit with spikes at 0.5 (bin 1), 1.5 and 1.7 (bin 2)
            ud = test_CellTypeClassifier_helpers.makeUnitData([0.5, 1.5, 1.7]);
            binned = CellTypeClassifier.binnedSpikeMatrix(ud, [0, 3], 1.0);
            tc.verifyEqual(binned(1, :), [1, 2, 0]);
        end

        function testBinnedMatrixNoSpikesIsZero(tc)
            ud = test_CellTypeClassifier_helpers.makeUnitData([]);
            binned = CellTypeClassifier.binnedSpikeMatrix(ud, [0, 5], 1.0);
            tc.verifyEqual(binned, zeros(1, 5));
        end

    end

    methods (Static, Access = private)

        function ud = makeUnitData(spike_times)
            ud = UnitData(1, 1, spike_times(:), zeros(82,1), "test_rec", 20000, 60);
        end

    end
end
