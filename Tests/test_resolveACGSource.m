classdef test_resolveACGSource < matlab.unittest.TestCase
% TEST_RESOLVEACGSOURCE  Harmonization.ACGSource decides which spikes enter the ACG.
%
% A unit fires regularly at 10 Hz in its baseline segment and in 5 ms doublets in a
% drug segment. "ACG" must use the baseline segment only (drug-free), "Concatenated"
% all segments.
%
% Run: runtests('test_resolveACGSource')

    methods (Test)

        function testACGUsesOwnSegmentOnly(tc)
            [ctc, ud, ph] = test_resolveACGSource.makeClassifier();
            ph.ACGSource = "ACG";
            acg = ctc.resolveACG(ud, ph);
            tc.verifyEqual(acg(test_resolveACGSource.binAt(ph, 0.005)), 0);
        end

        function testConcatenatedIncludesDrugSegment(tc)
            [ctc, ud, ph] = test_resolveACGSource.makeClassifier();
            ph.ACGSource = "Concatenated";
            acg = ctc.resolveACG(ud, ph);
            tc.verifyGreaterThan(acg(test_resolveACGSource.binAt(ph, 0.005)), 0);
        end

    end

    methods (Static)

        function [ctc, ud_base, ph] = makeClassifier()
            base_st = (0.05:0.1:60)';                    % no intervals < 100 ms
            dose_st = sort([base_st; base_st + 0.005]);  % doublets 5 ms apart
            ud_base = UnitData(1, 1, base_st, zeros(10, 1), "rec_base", 10000, 60);
            ud_dose = UnitData(1, 1, dose_st, zeros(10, 1), "rec_dose", 10000, 60);
            ud_base.UnitID = "u1";
            ud_dose.UnitID = "u1";
            ctc = CellTypeClassifier();
            ctc.Parameters = CellTypeClassifier.returnDefaultParams();
            ctc.OriginalUnitDataArray = [ud_base, ud_dose];
            ph = ctc.Parameters.Harmonization;
        end

        function b = binAt(ph, lag_s)
            % Index of the +lag_s bin in an ACG of 2*ACGLag/ACGBinSize + 1 bins
            b = round(ph.ACGLag / ph.ACGBinSize) + 1 + round(lag_s / ph.ACGBinSize);
        end

    end
end
