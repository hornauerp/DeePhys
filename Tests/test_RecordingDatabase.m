classdef test_RecordingDatabase < matlab.unittest.TestCase
% TEST_RECORDINGDATABASE  Registry round trip with a RecordingProcessor.
%
% Needs the Database Toolbox (MATLAB's sqlite interface); skipped otherwise.
%
% Run: runtests('test_RecordingDatabase')

    methods (TestClassSetup)
        function needDatabaseToolbox(tc)
            tc.assumeTrue(RecordingDatabase.isAvailable(), 'Database Toolbox not available.');
        end
    end

    methods (Test)

        function testRegisterProcessor(tc)
            sorting_path = make_mock_sorting(string(tempdir), 5, 50, 20000, 60);
            c1 = onCleanup(@() rmdir(sorting_path, 's'));
            sd = SpikeData.fromKilosort(sorting_path, struct('ChipID', 'DBChip', 'PlatingDate', '240101'));
            p.Analyses.Catch22 = false;
            p.QC.N_Units = 1;   % the mock sorting has only 5 templates
            proc = RecordingProcessor(sd, p);
            proc.runQC();

            f = string(tempname) + ".db";
            RecordingDatabase.resetInstance();
            db = RecordingDatabase.instance(f);
            c2 = onCleanup(@() test_RecordingDatabase.removeDB(f));
            tc.assertTrue(db.IsConnected);

            tc.verifyWarningFree(@() db.registerRecording(proc));
            tc.verifyWarningFree(@() db.registerUnits(proc));
            tc.verifyWarningFree(@() db.updateQCResults(proc));
            tc.verifyWarningFree(@() db.updateProcessingStatus(proc, 'runQC', 'completed', 1.5));

            T = db.queryRecordings(ChipID="DBChip");
            tc.verifyEqual(height(T), 1);
            tc.verifyEqual(double(T.n_good_units), numel(proc.Units));
            tc.verifyEqual(height(db.queryRecordings(Analysis="runQC")), 1);

            db.registerRecording(proc);   % re-registration updates, does not duplicate
            tc.verifyEqual(height(db.queryRecordings(ChipID="DBChip")), 1);
        end

    end

    methods (Static, Access = private)
        function removeDB(f)
            RecordingDatabase.resetInstance();
            if isfile(f), delete(f); end
        end
    end
end
