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
            [X_filt, ids_f] = fs.unitMatrix('all', string.empty, MinFiringRate=2.0);
            warning(warnState);
            tc.verifyLessThan(height(X_filt), height(X_all));
            % All remaining units should have FiringRate >= 2
            fr = fs.UnitTable.FiringRate(ismember(fs.UnitTable.UnitID, ids_f));
            tc.verifyTrue(all(fr >= 2.0));
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
