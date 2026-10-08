classdef test_computeConnectivityCCG < matlab.unittest.TestCase
% TEST_COMPUTECONNECTIVITYCCG  Unit tests for Functions/computeConnectivityCCG.m
%
% Run: runtests('test_computeConnectivityCCG')

    methods (Test)

        function testOutputStructFields(tc)
            [st, su] = tc.makeCorrelatedPair();
            params.BinSize = 0.001;
            params.Duration = 0.1;    % production default; 10-bin smoothing needs a CCG >= ~1.5x the kernel
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
            params.Duration = 0.1;    % production default; 10-bin smoothing needs a CCG >= ~1.5x the kernel
            params.Conv_w = 10;
            params.Alpha = 0.01;
            conn = computeConnectivityCCG(st, su, 3, 30000, params);
            tc.verifySize(conn.bd, [3, 3]);
        end

        function testCCGsHaveCorrectDimensions(tc)
            [st, su] = tc.makeCorrelatedPair();
            params.BinSize = 0.001;
            params.Duration = 0.1;    % production default; 10-bin smoothing needs a CCG >= ~1.5x the kernel
            params.Conv_w = 10;
            params.Alpha = 0.01;
            conn = computeConnectivityCCG(st, su, 3, 30000, params);
            tc.verifyEqual(size(conn.CCGs, 2), 3);
            tc.verifyEqual(size(conn.CCGs, 3), 3);
        end

        function testOneSidedSuppressionIsInhibitory(tc)
            % unit 2 is silenced 1-5 ms after each unit-1 spike: inhibitory 1 -> 2 only,
            % with either criterion
            [st, su] = tc.makeSuppressedPair(false);
            for crit = ["strict", "cellexplorer"]
                conn = tc.detect(st, su, crit);
                tc.verifyEqual(conn.InhibitoryConnection, [1 2], crit);
            end
        end

        function testSymmetricTroughIsNotInhibitory(tc)
            % two units that never fire within 5 ms of each other (e.g. a split unit):
            % the trough is symmetric, so neither direction is a synaptic inhibition
            [st, su] = tc.makeSuppressedPair(true);
            for crit = ["strict", "cellexplorer"]
                conn = tc.detect(st, su, crit);
                tc.verifyEmpty(conn.InhibitoryConnection, crit);
            end
        end

        function testReverseInhibitionOfExcitation(tc)
            % unit 2 drives unit 1 at 2 ms (excitatory 2 -> 1). The default "strict"
            % criterion reports no inhibition; the "cellexplorer" criterion also reports an
            % inhibitory 1 -> 2, the false positive CellExplorer warns about.
            rng(8);
            T = 300; st1 = sort(rand(3000, 1) * T); st2 = sort(rand(3000, 1) * T);
            k = rand(numel(st2), 1) < 0.3;
            st1 = sort([st1; st2(k) + 0.002 + 0.0003 * randn(sum(k), 1)]);
            st = [st1; st2]; su = [ones(numel(st1), 1); 2 * ones(numel(st2), 1)];
            [st, o] = sort(st); su = su(o);
            def = tc.detect(st, su);                                % function default
            tc.verifyTrue(ismember([2 1], def.ExcitatoryConnection, 'rows'));
            tc.verifyEmpty(def.InhibitoryConnection);
            ce = tc.detect(st, su, "cellexplorer");
            tc.verifyEqual(ce.ExcitatoryConnection, def.ExcitatoryConnection);
            tc.verifyTrue(ismember([1 2], ce.InhibitoryConnection, 'rows'));
        end

        function testProcessorDefaultIsStrict(tc)
            p = RecordingProcessor.returnDefaultParams();
            tc.verifyEqual(string(p.CCG.InhibitoryCriterion), "strict");
        end

    end

    methods (Static)
        function conn = detect(st, su, criterion)
            % criterion omitted: the field is not set, so the function default applies
            params = struct('BinSize', 0.001, 'Duration', 0.1, 'Conv_w', 10, 'Alpha', 0.001);
            if nargin >= 3, params.InhibitoryCriterion = criterion; end
            [~, conn] = evalc('computeConnectivityCCG(st, su, 2, 10000, params)');
        end

        function [spike_times, spike_units] = makeSuppressedPair(symmetric)
            rng(7);
            T = 300;
            st1 = sort(rand(3000, 1) * T);
            st2 = sort(rand(3000, 1) * T);
            d = st2 - st1';                                      % lag of unit 2 after unit 1
            if symmetric
                kill = any(abs(d) <= 0.005, 2);
            else
                kill = any(d > 0.001 & d <= 0.005, 2);
            end
            st2 = st2(~kill);
            spike_times = [st1; st2];
            spike_units = [ones(numel(st1), 1); 2 * ones(numel(st2), 1)];
            [spike_times, o] = sort(spike_times);
            spike_units = spike_units(o);
        end

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
