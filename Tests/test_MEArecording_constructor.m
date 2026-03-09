classdef test_MEArecording_constructor < matlab.unittest.TestCase
% TEST_MEARECORDING_CONSTRUCTOR  Smoke tests for MEArecording no-arg construction.
%
% These tests require only the class files to be on the MATLAB path.
% No recording data or external toolboxes are needed.
%
% Run: runtests('test_MEArecording_constructor')

    methods (Test)

        function testNoArgConstructorReturnsHandle(tc)
            obj = MEArecording();
            tc.verifyClass(obj, 'MEArecording');
            tc.verifyTrue(isa(obj, 'handle'));
        end

        function testNoArgConstructorPropertiesEmpty(tc)
            obj = MEArecording();
            tc.verifyEmpty(obj.Metadata);
            tc.verifyEmpty(obj.Units);
            tc.verifyEmpty(obj.RecordingInfo);
            tc.verifyEmpty(obj.Parameters);
        end

        function testEmptyStructEquivalentToNoArg(tc)
            obj1 = MEArecording();
            obj2 = MEArecording(struct(), struct());
            tc.verifyEqual(obj1.Metadata,     obj2.Metadata);
            tc.verifyEqual(obj1.Units,         obj2.Units);
            tc.verifyEqual(obj1.RecordingInfo, obj2.RecordingInfo);
        end

        function testArrayPreAllocation(tc)
            arr = MEArecording.empty(0, 3);
            tc.verifySize(arr, [0 3]);
            tc.verifyClass(arr, 'MEArecording');
        end

        function testFolderMethodDiscovery(tc)
            % Verifies that @folder method files are correctly resolved by MATLAB.
            % If a method body was accidentally left in the classdef AND in the
            % @folder file, MATLAB would use the classdef copy — this test
            % catches that by confirming method names are visible at all.
            m = methods('MEArecording');
            expected = {'inferConnectivityCCG', 'inferGraphFeatures', ...
                        'inferWaveformFeatures', 'detectBurstsISIN', ...
                        'generateUnits', 'getBurstStatistics'};
            for k = 1:numel(expected)
                tc.verifyTrue(ismember(expected{k}, m), ...
                    sprintf('"%s" must be discoverable as a method', expected{k}));
            end
        end

        function testIsHandleClass(tc)
            % MEArecording is a handle class — two variables pointing to the
            % same object should share state.
            obj1 = MEArecording();
            obj2 = obj1;
            obj1.NumUnitClusters = 42;
            tc.verifyEqual(obj2.NumUnitClusters, 42);
        end

    end
end
