classdef test_Unit_constructor < matlab.unittest.TestCase
% TEST_UNIT_CONSTRUCTOR  Smoke tests for Unit no-arg construction and static methods.
%
% No recording data or external toolboxes required.
%
% Run: runtests('test_Unit_constructor')

    methods (Test)

        function testNoArgConstructorReturnsHandle(tc)
            unit = Unit();
            tc.verifyClass(unit, 'Unit');
            tc.verifyTrue(isa(unit, 'handle'));
        end

        function testNoArgConstructorPropertiesEmpty(tc)
            unit = Unit();
            tc.verifyEmpty(unit.SpikeTimes);
            tc.verifyEmpty(unit.ReferenceWaveform);
            tc.verifyEmpty(unit.MEArecording);
            tc.verifyEmpty(unit.WaveformFeatures);
        end

        function testArrayPreAllocation(tc)
            arr = Unit.empty(0, 5);
            tc.verifySize(arr, [0 5]);
            tc.verifyClass(arr, 'Unit');
        end

        function testReturnFeatureNamesActivity(tc)
            names    = Unit.returnFeatureNames('act');
            expected = ["FiringRate", "MeanInterSpikeInterval", ...
                        "VarianceInterSpikeInterval", "CVInterSpikeInterval", ...
                        "PartialAutocorrelation"];
            tc.verifyEqual(names, expected);
        end

        function testReturnFeatureNamesRegularity(tc)
            names    = Unit.returnFeatureNames('reg');
            expected = ["RegularityFrequency", "RegularityMagnitude", "RegularityFit"];
            tc.verifyEqual(names, expected);
        end

        function testFolderMethodDiscovery(tc)
            m = methods('Unit');
            tc.verifyTrue(ismember('inferActivityFeatures', m), ...
                'inferActivityFeatures must be discoverable as a method');
            tc.verifyTrue(ismember('getRegularity', m), ...
                'getRegularity must be discoverable as a method');
        end

        function testIsHandleClass(tc)
            u1 = Unit();
            u2 = u1;
            u1.ClusterID = 99;
            tc.verifyEqual(u2.ClusterID, 99);
        end

        function testUnitIDDependentPropertyAccessible(tc)
            % unitID is a Dependent property — accessing it on an empty Unit
            % should not crash (it reads from MEArecording, which is empty here,
            % so the result may be empty but must not error).
            unit = Unit();
            tc.verifyWarningFree(@() unit.unitID);
        end

    end
end
