classdef test_computeConnectivityCCG < matlab.unittest.TestCase
% TEST_COMPUTECONNECTIVITYCCG  Unit tests for Functions/computeConnectivityCCG.m
%
% Run: runtests('test_computeConnectivityCCG')

    methods (Test)

        function testOutputStructFields(tc)
            [st, su] = tc.makeCorrelatedPair();
            params.BinSize = 0.001;
            params.Duration = 0.02;
            params.Conv_w = 10;
            params.Alpha = 0.01;
            conn = computeConnectivityCCG(st, su, 3, 30000, params);
            tc.verifyTrue(isfield(conn, 'ExcitatoryConnection'));
            tc.verifyTrue(isfield(conn, 'InhibitoryConnection'));
            tc.verifyTrue(isfield(conn, 'bd'));
            tc.verifyTrue(isfield(conn, 'CCGs'));
        end

        function testConnectivityMatrixSize(tc)
            [st, su] = tc.makeCorrelatedPair();
            params.BinSize = 0.001;
            params.Duration = 0.02;
            params.Conv_w = 10;
            params.Alpha = 0.01;
            conn = computeConnectivityCCG(st, su, 3, 30000, params);
            tc.verifySize(conn.bd, [3, 3]);
        end

        function testCCGsHaveCorrectDimensions(tc)
            [st, su] = tc.makeCorrelatedPair();
            params.BinSize = 0.001;
            params.Duration = 0.02;
            params.Conv_w = 10;
            params.Alpha = 0.01;
            conn = computeConnectivityCCG(st, su, 3, 30000, params);
            tc.verifyEqual(size(conn.CCGs, 2), 3);
            tc.verifyEqual(size(conn.CCGs, 3), 3);
        end

    end

    methods (Static)
        function [spike_times, spike_units] = makeCorrelatedPair()
            rng(42);
            n = 5000;
            st1 = sort(rand(n, 1) * 100);
            % Unit 2: copy of unit 1 with 2ms delay + jitter
            st2 = st1 + 0.002 + randn(n, 1) * 0.0005;
            st2 = st2(st2 > 0 & st2 < 100);
            % Unit 3: independent
            st3 = sort(rand(n, 1) * 100);
            spike_times = sort([st1; st2; st3]);
            spike_units = [ones(n, 1); 2*ones(length(st2), 1); 3*ones(n, 1)];
            [spike_times, idx] = sort(spike_times);
            spike_units = spike_units(idx);
        end
    end
end
