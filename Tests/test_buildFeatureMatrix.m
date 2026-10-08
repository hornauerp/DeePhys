classdef test_buildFeatureMatrix < matlab.unittest.TestCase
% TEST_BUILDFEATUREMATRIX  Waveform harmonisation across sampling rates.
%
% Run: runtests('test_buildFeatureMatrix')

    methods (Test)

        function testMixedSamplingRatesAlign(tc)
            % The same spike sampled at 10 and 20 kHz (the shorter one NaN-padded in a
            % common matrix) must give the same harmonised waveform when each unit
            % gets its own rate.
            ctc = CellTypeClassifier();
            ctc.Parameters = CellTypeClassifier.returnDefaultParams();
            spike = @(t) -exp(-(t / 0.25).^2) + 0.35 * exp(-((t - 0.6) / 0.35).^2);   % t in ms
            t10 = (-2:0.1:4)';    t20 = (-2:0.05:4)';
            wf  = nan(numel(t20), 2);
            wf(1:numel(t10), 1) = spike(t10);
            wf(:, 2)            = spike(t20);
            acg = rand(401, 2);
            [~, ~, A] = evalc_quiet(@() buildFeatureMatrix(ctc, wf, acg, [10000 20000]));
            tc.verifyLessThan(max(abs(A(:, 1) - A(:, 2))), 0.02);

            % A single rate for both (the former behaviour) stretches the 10 kHz unit
            [~, ~, B] = evalc_quiet(@() buildFeatureMatrix(ctc, wf, acg, 20000));
            tc.verifyGreaterThan(max(abs(B(:, 1) - B(:, 2))), 0.1);
        end

        function testBaselineACGNormalisation(tc)
            % "baseline": the outer fifth of the lags averages 1, independent of the rate
            ctc = CellTypeClassifier();
            ctc.Parameters = CellTypeClassifier.returnDefaultParams();
            ctc.Parameters.Harmonization.ACGNormalization = "baseline";
            lag = (-200:200)';
            shape = 1 - 0.9 * exp(-abs(lag) / 10);          % refractory dip, flat tail
            acg = [5 * shape, 50 * shape];                    % same shape, 10x the rate
            [~, ~, ~, A] = evalc_quiet(@() buildFeatureMatrix(ctc, randn(61, 2), acg, 10000));
            tail = [1:40, 362:401];
            tc.verifyEqual(mean(A(tail, :), 1), [1 1], 'AbsTol', 1e-12);
            tc.verifyEqual(A(:, 1), A(:, 2), 'AbsTol', 1e-12);
        end

        function testRateCountMustMatchUnits(tc)
            ctc = CellTypeClassifier();
            ctc.Parameters = CellTypeClassifier.returnDefaultParams();
            tc.verifyError(@() buildFeatureMatrix(ctc, randn(61, 3), rand(401, 3), [1e4 2e4]), ...
                'buildFeatureMatrix:samplingRate');
        end

    end
end

function varargout = evalc_quiet(fh)
% Call fh with its console output suppressed.
    [~, varargout{1:nargout}] = evalc('fh()');
end
