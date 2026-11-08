classdef test_computeWaveformFeatures < matlab.unittest.TestCase
% TEST_COMPUTEWAVEFORMFEATURES  Unit tests for Functions/computeWaveformFeatures.m
%
% Run: runtests('test_computeWaveformFeatures')

    methods (Test)

        function testOutputFieldNames(tc)
            wf = tc.makeSyntheticWaveform();
            feats = computeWaveformFeatures(wf, 10000);
            tc.verifyLength(feats, 1);
            expected = ["AUC_peak_1","AUC_trough","AUC_peak_2","Rise","Decay", ...
                        "HalfWidth","Asymmetry","T2Pdelay","T2Pratio"];
            actual = string(fieldnames(feats{1}))';
            tc.verifyEqual(actual, expected);
        end

        function testMultipleUnits(tc)
            wf1 = tc.makeSyntheticWaveform();
            wf2 = tc.makeSyntheticWaveform() * 0.8;
            feats = computeWaveformFeatures([wf1, wf2], 10000);
            tc.verifyLength(feats, 2);
        end

        function testHalfWidthPositive(tc)
            wf = tc.makeSyntheticWaveform();
            feats = computeWaveformFeatures(wf, 10000);
            tc.verifyGreaterThan(feats{1}.HalfWidth, 0);
        end

        function testT2PdelayPositive(tc)
            wf = tc.makeSyntheticWaveform();
            feats = computeWaveformFeatures(wf, 10000);
            tc.verifyGreaterThan(feats{1}.T2Pdelay, 0);
        end

        function testNoNaNsForGoodWaveform(tc)
            wf = tc.makeSyntheticWaveform();
            feats = computeWaveformFeatures(wf, 10000);
            vals = struct2array(feats{1});
            tc.verifyFalse(any(isnan(vals)));
        end

        function testRiseAndDecaySigns(tc)
            f = computeWaveformFeatures(tc.makeSyntheticWaveform(), 10000);
            tc.verifyGreaterThan(f{1}.Rise, 0);
            tc.verifyLessThan(f{1}.Decay, 0);
        end

        function testRiseAndDecayScaleWithAmplitude(tc)
            wf = tc.makeSyntheticWaveform();
            a = computeWaveformFeatures(wf, 10000);
            b = computeWaveformFeatures(0.5 * wf, 10000);
            tc.verifyEqual(b{1}.Rise, 0.5 * a{1}.Rise, 'RelTol', 1e-12);
            tc.verifyEqual(b{1}.Decay, 0.5 * a{1}.Decay, 'RelTol', 1e-12);
        end

        function testRiseAndDecayStableUnderRoundOff(tc)
            % Perturbations at round-off level change the slopes only at round-off level
            % (histogram state levels of slewrate could change them by tens of percent).
            rng(3);
            wf = tc.makeSyntheticWaveform();
            a = computeWaveformFeatures(wf, 10000);
            for k = 1:20
                b = computeWaveformFeatures(wf .* (1 + 1e-14 * randn(size(wf))), 10000);
                tc.verifyEqual(b{1}.Rise, a{1}.Rise, 'RelTol', 1e-9);
                tc.verifyEqual(b{1}.Decay, a{1}.Decay, 'RelTol', 1e-9);
            end
        end

        function testRiseOfALinearRamp(tc)
            % trough -1 at sample 21, linear to +0.5 at sample 31 (1 ms at 10 kHz), then flat:
            % the spline follows the ramp closely, so the 10-90 % slope is about 1.5 per ms
            down = linspace(0, -1, 11)'; up = linspace(-1, 0.5, 11)'; tail = linspace(0.5, -0.2, 20)';
            wf = [zeros(10, 1); down; up(2:end); 0.5 * ones(5, 1); tail(2:end)];
            f = computeWaveformFeatures(wf, 10000);
            tc.verifyEqual(f{1}.Rise, 1.5, 'RelTol', 0.05);
        end

    end

    methods (Static)
        function wf = makeSyntheticWaveform()
            t = linspace(0, 3, 60)';
            wf = 0.3 * exp(-((t - 1.0).^2) / 0.02) ...  % pre-peak
               - exp(-((t - 1.3).^2) / 0.03) ...          % trough
               + 0.4 * exp(-((t - 1.8).^2) / 0.05);      % post-peak
            wf = wf / max(abs(wf));
        end
    end
end
