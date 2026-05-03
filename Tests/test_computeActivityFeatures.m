classdef test_computeActivityFeatures < matlab.unittest.TestCase
% TEST_COMPUTEACTIVITYFEATURES  Unit tests for Functions/computeActivityFeatures.m
%
% Run: runtests('test_computeActivityFeatures')

    methods (Test)

        function testOutputColumnsMatch(tc)
            rng(42);
            spike_times = sort(rand(500, 1) * 100);
            params.FanoBinWidth = 0.1;
            params.RefractoryPeriod = 0.002;
            [act, ~, ~] = computeActivityFeatures(spike_times, 100, params);
            expected = ["FiringRate","MeanInterSpikeInterval","VarianceInterSpikeInterval", ...
                "CVInterSpikeInterval","CV2InterSpikeInterval","LocalVariation", ...
                "RevisedLocalVariation","FanoFactor","PartialAutocorrelation","ISIBimodality"];
            tc.verifyEqual(string(act.Properties.VariableNames), expected);
        end

        function testFiringRateCorrect(tc)
            spike_times = (1:100)';
            params.FanoBinWidth = 0.1;
            params.RefractoryPeriod = 0.002;
            [act, ~, ~] = computeActivityFeatures(spike_times, 100, params);
            tc.verifyEqual(act.FiringRate, 1.0, 'AbsTol', 1e-10);
        end

        function testInsufficientSpikesReturnsNaN(tc)
            spike_times = [1; 2];
            params.FanoBinWidth = 0.1;
            params.RefractoryPeriod = 0.002;
            [act, ~, ~] = computeActivityFeatures(spike_times, 100, params);
            tc.verifyTrue(all(isnan(act.Variables)));
        end

        function testSingleSpikeReturnsNaN(tc)
            spike_times = 5;
            [act, ~, ~] = computeActivityFeatures(spike_times, 100, struct());
            tc.verifyTrue(all(isnan(act.Variables)));
        end

        function testRegularityOptional(tc)
            rng(42);
            spike_times = sort(rand(200, 1) * 60);
            params.FanoBinWidth = 0.1;
            params.RefractoryPeriod = 0.002;
            [~, reg, ~] = computeActivityFeatures(spike_times, 60, params);
            tc.verifyTrue(isempty(reg) || width(reg) == 0);
        end

        function testRegularityReturnsTable(tc)
            rng(42);
            spike_times = sort(rand(200, 1) * 60);
            params.FanoBinWidth = 0.1;
            params.RefractoryPeriod = 0.002;
            params.Regularity.Binning = 0.1;
            params.Regularity.N_peaks = 3;
            [~, reg, ~] = computeActivityFeatures(spike_times, 60, params);
            tc.verifyEqual(width(reg), 3);
            expected = ["RegularityFrequency","RegularityMagnitude","RegularityFit"];
            tc.verifyEqual(string(reg.Properties.VariableNames), expected);
        end

        function testRegularSpikeTrainLowCV(tc)
            spike_times = (0.01:0.01:10)';
            params.FanoBinWidth = 0.1;
            params.RefractoryPeriod = 0.002;
            [act, ~, ~] = computeActivityFeatures(spike_times, 10, params);
            tc.verifyLessThan(act.CVInterSpikeInterval, 0.01);
        end

        function testOutputIsSingleRow(tc)
            rng(42);
            spike_times = sort(rand(300, 1) * 50);
            [act, ~, ~] = computeActivityFeatures(spike_times, 50, struct());
            tc.verifyEqual(height(act), 1);
        end

    end
end
