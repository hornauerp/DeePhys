classdef test_computeRegularity < matlab.unittest.TestCase
% TEST_COMPUTEREGULARITY  Unit tests for Functions/computeRegularity.m
%
% Run: runtests('test_computeRegularity')

    methods (Test)

        function testOutputSizesAreScalar(tc)
            binning = 0.01;
            t       = 0 : binning : 10 - binning;
            signal  = sin(2 * pi * 3 * t);
            [freq, mag, coeff] = computeRegularity(signal, binning, 3);
            tc.verifySize(freq,  [1 1]);
            tc.verifySize(mag,   [1 1]);
            tc.verifySize(coeff, [1 1]);
        end

        function testKnownFrequency1Hz(tc)
            % 10 s at 1 Hz tone, 10 ms bins → should detect ~1 Hz
            binning = 0.01;
            t       = 0 : binning : 10 - binning;
            signal  = sin(2 * pi * 1 * t);
            [frequency, magnitude, ~] = computeRegularity(signal, binning, 3);
            tc.verifyEqual(frequency, 1, 'AbsTol', 0.5);
            tc.verifyGreaterThan(magnitude, 0);
        end

        function testKnownFrequency5Hz(tc)
            binning = 0.01;
            t       = 0 : binning : 10 - binning;
            signal  = sin(2 * pi * 5 * t);
            [frequency, ~, ~] = computeRegularity(signal, binning, 3);
            tc.verifyEqual(frequency, 5, 'AbsTol', 0.5);
        end

        function testKnownFrequency10Hz(tc)
            binning = 0.01;
            t       = 0 : binning : 10 - binning;
            signal  = sin(2 * pi * 10 * t);
            [frequency, ~, ~] = computeRegularity(signal, binning, 3);
            tc.verifyEqual(frequency, 10, 'AbsTol', 1.0);
        end

        function testMagnitudeIsNonnegative(tc)
            binning = 0.01;
            t       = 0 : binning : 5 - binning;
            signal  = randn(size(t));
            [~, magnitude, ~] = computeRegularity(signal, binning, 3);
            tc.verifyGreaterThanOrEqual(magnitude, 0);
        end

        function testFitCoeffIsScalarForRealSignal(tc)
            % Well-behaved periodic signal should produce a finite fit_coeff
            binning = 0.01;
            t       = 0 : binning : 20 - binning;
            signal  = abs(sin(2 * pi * 2 * t)) + 0.1 * randn(size(t));
            [~, ~, fit_coeff] = computeRegularity(signal, binning, 5);
            tc.verifySize(fit_coeff, [1 1]);
        end

        function testFlatSignalProducesScalarCoeff(tc)
            % Flat signal → FFT all-zero after DC removal → fit may fail → NaN
            % Key requirement: output is still a scalar (not an error)
            signal  = ones(1, 500);
            binning = 0.01;
            [~, ~, fit_coeff] = computeRegularity(signal, binning, 3);
            tc.verifySize(fit_coeff, [1 1]);
        end

        function testFrequencyNonnegative(tc)
            binning = 0.005;
            t       = 0 : binning : 8 - binning;
            signal  = randn(size(t));
            [frequency, ~, ~] = computeRegularity(signal, binning, 4);
            tc.verifyGreaterThanOrEqual(frequency, 0);
        end

    end
end
