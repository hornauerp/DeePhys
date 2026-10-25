classdef test_exportReference < matlab.unittest.TestCase
% TEST_EXPORTREFERENCE  Golden-file exporter (Tests/export_reference.m): files, metadata, 0-based indices.
%
% Run: runtests('test_exportReference')

    properties
        OutDir
    end

    methods (TestClassSetup)
        function exportOnce(tc)
            tc.OutDir = fullfile(tempname, 'golden');
            tc.addTeardown(@() rmdir(fileparts(tc.OutDir), 's'));
            evalc('export_reference(tc.OutDir)');
        end
    end

    methods (Test)

        function testSourceAndReadmeWritten(tc)
            src = jsondecode(fileread(fullfile(tc.OutDir, 'SOURCE.json')));
            tc.verifyNotEmpty(src.deephys_commit);
            tc.verifyNotEmpty(src.matlab_version);
            tc.verifyGreaterThan(src.n_cases, 50);
            tc.verifyEqual(numel(src.cases), src.n_cases);
            tc.verifyTrue(isfile(fullfile(tc.OutDir, 'README.md')));
        end

        function testExpectedCaseFilesExist(tc)
            for rel = ["CCG/trio_bin1ms_dur100ms_counts.h5", "sttc/corr_dt50ms.h5", ...
                    "detectNetworkBursts/network.h5", "buildKNNGraph/three_clusters_k10.h5", ...
                    "ComBat/two_batches.h5", "clustering_coef_bd/random_n24.h5"]
                tc.verifyTrue(isfile(fullfile(tc.OutDir, rel)), rel);
            end
            tc.verifyTrue(isfile(fullfile(tc.OutDir, 'communityEnrichmentTest', 'pooled_T.parquet')));
        end

        function testFileMetadata(tc)
            f = fullfile(tc.OutDir, 'CCG', 'trio_bin1ms_dur100ms_counts.h5');
            tc.verifyEqual(string(h5readatt(f, '/', 'layout')), "matlab-column-major");
            opts = jsondecode(h5readatt(f, '/', 'options'));
            tc.verifyEqual(opts.binSize, 0.001);
            ccg = h5read(f, '/outputs/ccg');
            tc.verifySize(ccg, [101 3 3]);
        end

        function testIndicesAreZeroBased(tc)
            f = fullfile(tc.OutDir, 'computeConnectivityCCG', 'trio_strict.h5');
            su = h5read(f, '/inputs/spike_units');
            tc.verifyEqual(min(su), 0);
            tc.verifyEqual(max(su), 2);
            exc = h5read(f, '/outputs/ExcitatoryConnection');
            tc.verifyEqual(exc, [0 1]);           % unit 1 drives unit 2 in MATLAB numbering
            idx = jsondecode(h5readatt(f, '/', 'index_valued'));
            tc.verifyTrue(ismember('inputs/spike_units', idx));
            seeds = h5read(fullfile(tc.OutDir, 'propagateLabels', 'balanced_seeds.h5'), '/inputs/seed_idx');
            tc.verifyEqual(min(seeds(:)), 0);
        end

        function testNaNSurvives(tc)
            v = h5read(fullfile(tc.OutDir, 'sttc', 'single_spike_vs_two.h5'), '/outputs/tileCoef');
            tc.verifyTrue(isscalar(v));
            nanf = fullfile(tc.OutDir, 'computeRegularity', 'poisson_train_bin100ms_peaks3.h5');
            tc.verifyTrue(isnan(h5read(nanf, '/outputs/fit_coeff')));
        end

    end
end
