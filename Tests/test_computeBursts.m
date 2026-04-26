classdef test_computeBursts < matlab.unittest.TestCase
% TEST_COMPUTEBURSTS  Unit tests for Functions/computeBursts.m
%
% Run: runtests('test_computeBursts')

    methods (Test)

        function testLogISIDetectsBursts(tc)
            [st, su] = tc.makeBurstyData();
            n_units = max(su);
            fr  = ones(n_units, 1) * 5;
            cv  = ones(n_units, 1) * 1.2;
            params.Method = "logISI";
            params.Outlier = struct('Method', '', 'ThresholdFactor', 3);
            [bursts, features] = computeBursts(st, su, fr, cv, n_units, max(st) + 1, params);
            tc.verifyFalse(isempty(bursts.T_start));
            tc.verifyClass(features, 'table');
        end

        function testBurstFeaturesTableColumns(tc)
            [st, su] = tc.makeBurstyData();
            n_units = max(su);
            fr  = ones(n_units, 1) * 5;
            cv  = ones(n_units, 1) * 1.2;
            params.Method = "logISI";
            params.Outlier = struct('Method', '', 'ThresholdFactor', 3);
            [~, features] = computeBursts(st, su, fr, cv, n_units, max(st) + 1, params);
            expected = ["MeanInterBurstInterval","BurstDuration","RiseTime","FallTime", ...
                "PeakFR","StdFR","StdBurstDuration","StdInterBurstInterval", ...
                "IntraFiringRate","InterFiringRate"];
            tc.verifyEqual(string(features.Properties.VariableNames), expected);
        end

        function testNoBurstsInRandomSpikes(tc)
            rng(42);
            st = sort(rand(50, 1) * 100);
            su = randi(3, 50, 1);
            fr  = [2; 2; 2];
            cv  = [1.5; 1.5; 1.5];
            params.Method = "logISI";
            params.Outlier = struct('Method', '', 'ThresholdFactor', 3);
            [bursts, features] = computeBursts(st, su, fr, cv, 3, 100, params);
            tc.verifyTrue(isempty(bursts.T_start) || all(isnan(features.Variables)));
        end

        function testOutputStructFields(tc)
            [st, su] = tc.makeBurstyData();
            n_units = max(su);
            fr  = ones(n_units, 1) * 5;
            cv  = ones(n_units, 1) * 1.2;
            params.Method = "logISI";
            params.Outlier = struct('Method', '', 'ThresholdFactor', 3);
            [bursts, ~] = computeBursts(st, su, fr, cv, n_units, max(st) + 1, params);
            tc.verifyTrue(isfield(bursts, 'T_start'));
            tc.verifyTrue(isfield(bursts, 'T_end'));
            tc.verifyTrue(isfield(bursts, 'Method'));
        end

    end

    methods (Static)
        function [spike_times, spike_units] = makeBurstyData()
            rng(42);
            burst_centers = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95];
            st = [];
            su = [];
            for b = 1:length(burst_centers)
                n_spikes = randi([30, 60]);
                burst_spikes = burst_centers(b) + randn(n_spikes, 1) * 0.05;
                st = [st; burst_spikes]; %#ok<AGROW>
                su = [su; randi(5, n_spikes, 1)]; %#ok<AGROW>
            end
            % Add sparse inter-burst spikes
            inter_spikes = sort(rand(50, 1) * 100);
            st = [st; inter_spikes];
            su = [su; randi(5, 50, 1)];
            [st, idx] = sort(st);
            su = su(idx);
        end
    end
end
