classdef test_parseStructParameters < matlab.unittest.TestCase
% TEST_PARSESTRUCTPARAMETERS  Unit tests for Functions/parseStructParameters.m
%
% Run: runtests('test_parseStructParameters')

    methods (Test)

        function testEmptyOverridesReturnsDefaultsUnchanged(tc)
            defaults.A.x = 1;
            defaults.A.y = 2;
            overrides    = struct();
            merged       = parseStructParameters(defaults, overrides);
            tc.verifyEqual(merged, defaults);
        end

        function testValidOverrideUpdatesSubfield(tc)
            defaults.A.x = 1;
            defaults.A.y = 2;
            overrides.A.x = 99;
            merged = parseStructParameters(defaults, overrides);
            tc.verifyEqual(merged.A.x, 99);
            tc.verifyEqual(merged.A.y, 2);   % unchanged
        end

        function testMultipleTopLevelFieldsUpdated(tc)
            defaults.A.x = 1;
            defaults.A.y = 2;
            defaults.B.z = 3;
            overrides.A.x = 10;
            overrides.B.z = 30;
            merged = parseStructParameters(defaults, overrides);
            tc.verifyEqual(merged.A.x, 10);
            tc.verifyEqual(merged.B.z, 30);
            tc.verifyEqual(merged.A.y, 2);   % unchanged
        end

        function testUnknownTopLevelFieldEmitsWarning(tc)
            defaults.A.x = 1;
            overrides.UnknownTop.x = 99;
            lastwarn('', '');
            parseStructParameters(defaults, overrides);
            [msg, ~] = lastwarn();
            tc.verifyNotEmpty(msg, ...
                'Expected a warning for unrecognised top-level field');
        end

        function testUnknownSubfieldEmitsWarning(tc)
            defaults.A.x = 1;
            overrides.A.unknown_sub = 99;
            lastwarn('', '');
            parseStructParameters(defaults, overrides);
            [msg, ~] = lastwarn();
            tc.verifyNotEmpty(msg, ...
                'Expected a warning for unrecognised sub-level field');
        end

        function testNoWarningForValidOverride(tc)
            defaults.A.x = 1;
            overrides.A.x = 2;
            lastwarn('', '');
            parseStructParameters(defaults, overrides);
            [msg, ~] = lastwarn();
            tc.verifyEmpty(msg, ...
                'No warning expected when all override fields are recognised');
        end

        function testDefaultsStructNotMutated(tc)
            % MATLAB passes structs by value — original defaults must be unchanged
            defaults.A.x = 1;
            overrides.A.x = 99;
            parseStructParameters(defaults, overrides);
            tc.verifyEqual(defaults.A.x, 1);
        end

        function testStringValueOverride(tc)
            defaults.Outputs.Format = "csv";
            overrides.Outputs.Format = "mat";
            merged = parseStructParameters(defaults, overrides);
            tc.verifyEqual(merged.Outputs.Format, "mat");
        end

        function testLogicalValueOverride(tc)
            defaults.Analyses.ComputeBursts = false;
            overrides.Analyses.ComputeBursts = true;
            merged = parseStructParameters(defaults, overrides);
            tc.verifyTrue(merged.Analyses.ComputeBursts);
        end

    end
end
