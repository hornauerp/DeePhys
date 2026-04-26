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
