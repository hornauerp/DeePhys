classdef test_parseDateRobust < matlab.unittest.TestCase
% TEST_PARSEDATEROBUST  Date formats used in metadata and documentation.
%
% Run: runtests('test_parseDateRobust')

    methods (Test)

        function testCompactFormats(tc)
            tc.verifyEqual(parseDateRobust('221207'), datetime(2022, 12, 7));
            tc.verifyEqual(parseDateRobust(221207), datetime(2022, 12, 7));
            tc.verifyEqual(parseDateRobust("20221207"), datetime(2022, 12, 7));
            tc.verifyEqual(parseDateRobust('991231'), datetime(1999, 12, 31));   % MATLAB pivot
        end

        function testIsoAndSeparators(tc)
            tc.verifyEqual(parseDateRobust('2024-01-15'), datetime(2024, 1, 15));
            tc.verifyEqual(parseDateRobust('2024/1/5'), datetime(2024, 1, 5));
            tc.verifyEqual(parseDateRobust(datetime(2024, 3, 2)), datetime(2024, 3, 2));
        end

        function testInvalidGivesNaTWithWarning(tc)
            dt = tc.verifyWarning(@() parseDateRobust('2024-13-01'), 'DeePhys:UnknownDateFormat');
            tc.verifyTrue(isnat(dt));
            dt = tc.verifyWarning(@() parseDateRobust('15.01.24'), 'DeePhys:UnknownDateFormat');
            tc.verifyTrue(isnat(dt));
        end

        function testDIVFromIsoDates(tc)
            tc.verifyEqual(days(parseDateRobust('2024-01-29') - parseDateRobust('240101')), 28);
        end

    end
end
