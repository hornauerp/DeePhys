classdef test_calibration < matlab.unittest.TestCase
% TEST_CALIBRATION  Control-calibrated responder calls and label-free ACG smoothing.
%
% Run: runtests('test_calibration')

    methods (Test)

        function testResponseTableAndCalibratedCalls(tc)
            [ud, meta] = test_calibration.makeExperiment();
            p = struct('BinSize', 20, 'HighDoses', 2);
            T = CellTypeClassifier.responseTable(ud, meta, ["ChipID", "PlatingDate"], ...
                "Concentration", 0, ["A_1"; "B_1"], p);
            tc.verifyEqual(height(T), 100);                        % 50 tracked units per culture
            resp = startsWith(T.UnitID, "A") & ismember(extractAfter(T.UnitID, "_u"), string(1:5));
            tc.verifyEqual(median(T.LFC(resp)), 2, 'AbsTol', 0.3);  % 4x rate -> log2 = 2
            tc.verifyLessThan(abs(median(T.LFCc(~resp))), 0.1);
            is_ctrl = startsWith(T.UnitID, "B");
            [inc, ~, info] = CellTypeClassifier.calibrateResponses(T, is_ctrl, 0.05, 0.01);
            tc.verifyTrue(all(inc(resp)));
            tc.verifyEqual(info.thr_inc, quantile(T.LFCc(is_ctrl), 0.95), 'AbsTol', 1e-12);
            tc.verifyLessThan(info.fdr_inc, 0.5);
        end

        function testTargetFDRThreshold(tc)
            % treated: 900 null units like the controls + 100 responders with LFCc ~ 3
            rng(5);
            ctrl = 0.5 * randn(1000, 1);
            trt  = [0.5 * randn(900, 1); 3 + 0.5 * randn(100, 1)];
            LFCc = [trt; ctrl];
            T = table(LFCc, zeros(2000, 1), 'VariableNames', {'LFCc', 'p'});
            is_ctrl = [false(1000, 1); true(1000, 1)];
            [inc_strict, ~, i1] = CellTypeClassifier.calibrateResponses(T, is_ctrl, 0.001, 0.01);
            [inc_fdr,    ~, i2] = CellTypeClassifier.calibrateResponses(T, is_ctrl, 0.001, 0.01, 0.2);
            tc.verifyLessThanOrEqual(i2.fdr_inc, 0.2 + 1e-12);
            tc.verifyGreaterThanOrEqual(sum(inc_fdr & ~is_ctrl), sum(inc_strict & ~is_ctrl));
            tc.verifyLessThan(i2.thr_inc, i1.thr_inc);
        end

        function testExternalCutoffShiftedScores(tc)
            % excitatory units shifted above 0.5 (as for external units attached via interneurons)
            rng(6);
            p_exc = 1 ./ (1 + exp(-(0.2 + 0.4 * randn(60, 1))));    % ~0.55
            p_inh = 1 ./ (1 + exp(-(1.9 + 0.15 * randn(340, 1))));  % ~0.87
            [t, info] = CellTypeClassifier.externalCutoff([p_exc; p_inh]);
            tc.verifyTrue(info.bimodal);
            tc.verifyGreaterThan(t, max(prctile(p_exc, 90), 0.5));
            tc.verifyLessThan(t, prctile(p_inh, 10));
            tc.verifyEqual(info.inhibitory_share, 340 / 400, 'AbsTol', 0.05);
        end

        function testExternalCutoffUnimodalKeepsHalf(tc)
            rng(8);
            p = 1 ./ (1 + exp(-0.3 * randn(300, 1)));
            [t, info] = CellTypeClassifier.externalCutoff(p);
            tc.verifyEqual(t, 0.5);
            tc.verifyFalse(info.bimodal);
            [t2, info2] = CellTypeClassifier.externalCutoff(rand(10, 1));
            tc.verifyEqual(t2, 0.5);
            tc.verifyEqual(info2.method, "fixed (too few units)");
        end

        function testCalibrationNeedsControls(tc)
            T = table(["a"; "b"], ["c"; "c"], ["r"; "r"], [1; 1], [1; 1], [0; 1], [0; 1], [0.5; 0.001], ...
                'VariableNames', {'UnitID', 'Culture', 'BaselineRecordingID', 'FR0', 'FRhigh', 'LFC', 'LFCc', 'p'});
            tc.verifyError(@() CellTypeClassifier.calibrateResponses(T, false(2, 1), 0.01, 0.01), ...
                'CellTypeClassifier:noControls');
        end

        function testSmoothingKeepsCentreAndMass(tc)
            acg = zeros(201, 2); acg(101, :) = 5; acg(120, 1) = 10; acg(60:140, 2) = 1;
            s = CellTypeClassifier.smoothACGs(acg, 1, 0.0005);
            tc.verifyEqual(s(101, :), acg(101, :));
            tc.verifyEqual(sum(s([1:100, 102:end], 1)), 10, 'RelTol', 0.02);
            tc.verifyGreaterThan(s(118, 1), 0);
            tc.verifyEqual(CellTypeClassifier.smoothACGs(acg, 0, 0.0005), acg);
        end

        function testSmoothingWidthShrinksWithData(tc)
            % sparse trains need more smoothing than dense ones (same underlying process)
            rng(4);
            mk = @(rate, T) cumsum(0.002 + exprnd(1 / rate, round(rate * T), 1));
            sparse_tr = arrayfun(@(i) mk(2, 600), 1:60, 'UniformOutput', false);
            dense_tr  = arrayfun(@(i) mk(20, 600), 1:60, 'UniformOutput', false);
            s_sparse = CellTypeClassifier.acgSmoothingWidth(sparse_tr, [], 0.0005, 0.05);
            s_dense  = CellTypeClassifier.acgSmoothingWidth(dense_tr, [], 0.0005, 0.05);
            tc.verifyGreaterThan(s_sparse, s_dense);
        end

    end

    methods (Static)

        function [ud, meta] = makeExperiment()
            % Two cultures (A treated, B control), recordings at 0, 1, 10 (dose), 600 s each.
            % 50 units per culture at 2-5 Hz; units 1-5 of culture A fire 4x faster at doses > 0.
            rng(7);
            chips = ["A", "B"]; conc = [0 1 10]; dur = 600;
            RecordingID = strings(0, 1); ChipID = RecordingID; PlatingDate = RecordingID; Concentration = [];
            ud = UnitData.empty(1, 0);
            for c = chips
                base_rate = 2 + 3 * rand(50, 1);
                for d = conc
                    rid = c + "_" + d;
                    RecordingID(end+1, 1) = rid; ChipID(end+1, 1) = c; PlatingDate(end+1, 1) = "1"; %#ok<AGROW>
                    Concentration(end+1, 1) = d; %#ok<AGROW>
                    for u = 1:50
                        r = base_rate(u);
                        if c == "A" && u <= 5 && d > 0, r = 4 * r; end
                        st = sort(rand(poissrnd(r * dur), 1) * dur);
                        x = UnitData(u, u, st, zeros(10, 1), rid, 10000, dur);
                        x.UnitID = c + "_u" + u;
                        ud(end+1) = x; %#ok<AGROW>
                    end
                end
            end
            meta = table(RecordingID, ChipID, PlatingDate, Concentration);
        end

    end
end
