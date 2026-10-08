classdef test_FeatureStore < matlab.unittest.TestCase
% TEST_FEATURESTORE  Unit tests for FeatureStore table assembly and subsetting.
%
% Run: runtests('test_FeatureStore')

    methods (Test)

        function testSubsetByMetadataFiltersCorrectly(tc)
            fs = test_FeatureStore.makeMinimalFS();
            fs2 = fs.subsetByMetadata('Genotype', {'WT'});
            tc.verifyTrue(all(string(fs2.MetadataTable.Genotype) == "WT"));
            tc.verifyTrue(all(string(fs2.UnitTable.Genotype) == "WT"));
        end

        function testSubsetByMetadataNumeric(tc)
            fs = test_FeatureStore.makeMinimalFS();
            fs2 = fs.subsetByMetadata('DIV', 14);
            tc.verifyEqual(height(fs2.MetadataTable), 1);
        end

        function testSubsetByMetadataMissingFieldErrors(tc)
            fs = test_FeatureStore.makeMinimalFS();
            tc.verifyError( ...
                @() fs.subsetByMetadata('NonExistentField', 'x'), ...
                'FeatureStore:subsetByMetadata');
        end

        function testSubsetByRecordingIDsFiltersTables(tc)
            fs = test_FeatureStore.makeMinimalFS();
            ids = unique(fs.MetadataTable.RecordingID);
            fs2 = fs.subset(ids(1));
            tc.verifyEqual(height(fs2.MetadataTable), 1);
            tc.verifyTrue(all(fs2.UnitTable.RecordingID == ids(1)));
        end

        function testUnitMatrixExcludesIDColumns(tc)
            fs = test_FeatureStore.makeMinimalFS();
            [X, unit_ids] = fs.unitMatrix('all');
            % Identity columns must always be excluded from the feature matrix
            col_names = string(X.Properties.VariableNames);
            tc.verifyFalse(ismember('RecordingID', col_names));
            tc.verifyFalse(ismember('UnitID', col_names));
            % Feature columns must be present
            tc.verifyTrue(ismember('FiringRate', col_names));
            tc.verifyEqual(numel(unit_ids), height(X));
            tc.verifyEqual(height(X), height(fs.UnitTable));
        end

        function testUnitMatrixMinFiringRateFilter(tc)
            fs = test_FeatureStore.makeMinimalFS();
            [X_all, ~] = fs.unitMatrix('all');
            warnState  = warning('off', 'FeatureStore:unitMatrix');
            [X_filt, ids_f, row_idx] = fs.unitMatrix('all', string.empty, MinFiringRate=2.0);
            warning(warnState);
            tc.verifyLessThan(height(X_filt), height(X_all));
            % All remaining units should have FiringRate >= 2, and row_idx must align
            % X rows with UnitTable rows (UnitIDs repeat across recordings).
            tc.verifyEqual(numel(row_idx), height(X_filt));
            tc.verifyEqual(fs.UnitTable.UnitID(row_idx), ids_f);
            tc.verifyTrue(all(fs.UnitTable.FiringRate(row_idx) >= 2.0));
        end

        function testRecordingMatrixRowsMatchMetadata(tc)
            fs = test_FeatureStore.makeMinimalFS();
            [X, rec_ids] = fs.recordingMatrix('all');
            tc.verifyEqual(numel(rec_ids), height(fs.MetadataTable));
            tc.verifyEqual(height(X), numel(rec_ids));
        end

        function testSaveLoadRoundTrip(tc)
            fs = test_FeatureStore.makeMinimalFS();
            tmp = [tempname '.mat'];
            cleanup = onCleanup(@() delete(tmp));
            fs.save(tmp);
            fs2 = FeatureStore.load(tmp);
            tc.verifyEqual(height(fs2.UnitTable),     height(fs.UnitTable));
            tc.verifyEqual(height(fs2.RecordingTable), height(fs.RecordingTable));
            tc.verifyEqual(height(fs2.MetadataTable),  height(fs.MetadataTable));
        end

        function testBuildCultureIDs(tc)
            meta = table( ...
                ["r1";"r2";"r3"], ...
                ["A";"A";"B"], ...
                [1;1;2], ...
                'VariableNames', {'RecordingID','ChipID','PlatingDate'});
            cids = FeatureStore.buildCultureIDs(meta, ["ChipID","PlatingDate"]);
            tc.verifyEqual(cids(1), cids(2));  % same ChipID+PlatingDate
            tc.verifyNotEqual(cids(1), cids(3));
        end

        function testNetworkSelectorsMatchStoredNames(tc)
            % Network regularity is stored as NetworkRegularity*, network catch22 under the
            % plain catch22 names (SC_ marks single-unit catch22).
            fs = test_FeatureStore.makeMinimalFS();
            fs.RecordingTable.NetworkRegularityFrequency = [1; 2];
            fs.RecordingTable.NetworkRegularityMagnitude = [3; 4];
            fs.RecordingTable.NetworkRegularityFit       = [5; 6];
            fs.RecordingTable.DN_HistogramMode_5         = [7; 8];
            fs.RecordingTable.CO_f1ecac                  = [9; 10];
            X = fs.recordingMatrix('NetworkRegularity');
            tc.verifyEqual(size(X, 2), 3);
            X = fs.recordingMatrix('NetworkCatch22');
            tc.verifyEqual(size(X, 2), 2);
        end

        function testACGColumnsAndParentSubstitution(tc)
            fs = test_FeatureStore.makeMinimalFS();
            n = height(fs.UnitTable);
            fs.UnitTable.ACG1 = ones(n, 1);         fs.UnitTable.ACG2 = 2 * ones(n, 1);
            fs.UnitTable.Parent_ACG1 = 10 * ones(n, 1); fs.UnitTable.Parent_ACG2 = 20 * ones(n, 1);
            X_child = fs.unitMatrix('ACG');
            tc.verifyEqual(X_child{:, :}, repmat([1 2], n, 1));
            X_parent = fs.unitMatrix('ACG', 'ACG');
            tc.verifyEqual(X_parent{:, :}, repmat([10 20], n, 1));
            tc.verifyEqual(string(X_parent.Properties.VariableNames), ["ACG1", "ACG2"]);
        end

        function testWithoutRecordingFeatures(tc)
            % Processors without network features give an empty RecordingTable (no columns).
            fs = test_FeatureStore.makeMinimalFS();
            fs.RecordingTable = table();
            fs2 = fs.subset("r1");
            tc.verifyEqual(height(fs2.UnitTable), 3);
            tc.verifyError(@() fs.recordingMatrix('all'), 'FeatureStore:noRecordingFeatures');
        end

        function testSubsetOnEmptyResultIsHandled(tc)
            fs = test_FeatureStore.makeMinimalFS();
            % Subsetting by a non-existent recording ID should return empty tables
            fs2 = fs.subset("nonexistent_id");
            tc.verifyEqual(height(fs2.UnitTable), 0);
            tc.verifyEqual(height(fs2.MetadataTable), 0);
        end

    end

    methods (Static, Access = private)

        function fs = makeMinimalFS()
            % Build a minimal FeatureStore from two synthetic recordings.
            unit_table = table( ...
                ["u1";"u2";"u3";"u4";"u5";"u6"], ...
                ["r1";"r1";"r1";"r2";"r2";"r2"], ...
                ["WT";"WT";"WT";"KO";"KO";"KO"], ...
                [14;14;14;21;21;21], ...
                [1.5;3.0;0.5;2.0;4.0;1.0], ...
                [0.8;0.4;1.2;0.6;0.3;0.9], ...
                'VariableNames', {'UnitID','RecordingID','Genotype','DIV','FiringRate','CVInterSpikeInterval'});

            rec_table = table( ...
                ["r1";"r2"], ...
                ["WT";"KO"], ...
                [14;21], ...
                [2.0;1.5], ...
                'VariableNames', {'RecordingID','Genotype','DIV','BurstRate'});

            meta_table = table( ...
                ["r1";"r2"], ...
                ["WT";"KO"], ...
                [14;21], ...
                'VariableNames', {'RecordingID','Genotype','DIV'});

            fs = FeatureStore();
            fs.UnitTable      = unit_table;
            fs.RecordingTable = rec_table;
            fs.MetadataTable  = meta_table;
        end

    end
end
