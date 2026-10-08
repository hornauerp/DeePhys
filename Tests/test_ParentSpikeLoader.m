classdef test_ParentSpikeLoader < matlab.unittest.TestCase
% TEST_PARENTSPIKELOADER  Parent ACGs are keyed by the right TemplateID for 0- and 1-based exports.
%
% Template 1 fires regularly at 10 Hz (no intervals < 100 ms), template 2 in 5 ms doublets.
% TemplateID 2 (SpikeData convention: Kilosort value + 1) must get the doublet ACG whether
% the parent folder stores 0-based (Kilosort) or 1-based (qc_output export) values.
%
% Run: runtests('test_ParentSpikeLoader')

    methods (Test)

        function testZeroBasedKilosort(tc)
            acg_map = test_ParentSpikeLoader.load(tc, 0);
            test_ParentSpikeLoader.verifyTemplates(tc, acg_map);
        end

        function testOneBasedExport(tc)
            acg_map = test_ParentSpikeLoader.load(tc, 1);
            test_ParentSpikeLoader.verifyTemplates(tc, acg_map);
        end

    end

    methods (Static)

        function acg_map = load(tc, base)
            d = tempname; mkdir(d);
            tc.addTeardown(@() rmdir(d, 's'));
            sr = 10000;
            regular  = (0.05:0.1:60)';
            doublets = sort([(0.07:0.2:60)'; (0.07:0.2:60)' + 0.005]);
            noise    = sort(rand(300, 1) * 60);
            t   = [regular; doublets; noise];
            ids = [zeros(size(regular)); ones(size(doublets)); 2 * ones(size(noise))] + base;
            [t, o] = sort(t);
            writeNPY(round(t * sr), fullfile(d, 'spike_times.npy'));
            writeNPY(ids(o), fullfile(d, 'spike_templates.npy'));
            writeNPY(zeros(3, 4, 2), fullfile(d, 'templates.npy'));
            ParentSpikeLoader.clearCache();
            acg_map = ParentSpikeLoader.loadFullACGs(string(d), sr, struct('BinSize', 0.0005, 'Lag', 0.1));
        end

        function verifyTemplates(tc, acg_map)
            tc.verifyEqual(sort(string(keys(acg_map))), ["1", "2", "3"]);
            b5 = round(0.1 / 0.0005) + 1 + round(0.005 / 0.0005);   % +5 ms bin
            a1 = acg_map('1'); a2 = acg_map('2');
            tc.verifyEqual(a1(b5), 0);           % regular train: no 5 ms intervals
            tc.verifyGreaterThan(a2(b5), 0.5);   % doublets: peak at 5 ms
        end

    end
end
