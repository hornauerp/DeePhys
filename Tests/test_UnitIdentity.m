classdef test_UnitIdentity < matlab.unittest.TestCase
% TEST_UNITIDENTITY  UnitIDs: shared by the segments of one sorting, distinct for separately
% sorted recordings of the same culture, independent of the mount point.
%
% Run: runtests('test_UnitIdentity')

    properties
        Root
    end

    methods (TestMethodSetup)
        function makeRoot(tc)
            tc.Root = string(tempname);
            mkdir(tc.Root);
            tc.addTeardown(@() rmdir(tc.Root, 's'));
        end
    end

    methods (Test)

        function testSegmentsOfOneSortingShareUnitIDs(tc)
            ks = make_mock_sorting(tc.Root, 5, 80, 20000, 60);
            s0 = makeSegment(ks, 'segment_0', 1);
            s1 = makeSegment(ks, 'segment_1', 2);
            [p0, sd0] = processed(s0);
            [p1, sd1] = processed(s1);
            tc.verifyEqual(sd0.SortingID, sd1.SortingID);
            tc.verifyNotEqual(sd0.RecordingID, sd1.RecordingID);
            [common, i0, i1] = intersect([p0.Units.TemplateID], [p1.Units.TemplateID]);
            tc.verifyNotEmpty(common);
            tc.verifyEqual(string({p0.Units(i0).UnitID}), string({p1.Units(i1).UnitID}));
        end

        function testSeparateSortingsOfOneCultureDiffer(tc)
            % same chip, same templates and electrodes, but another day's sorting
            ks = make_mock_sorting(tc.Root, 5, 80, 20000, 60);
            other = fullfile(tc.Root, "day2", "sorter_output");
            mkdir(fileparts(other)); copyfile(ks, other);
            t = readNPY(fullfile(other, 'spike_times.npy')); t(end) = t(end) + 1;
            writeNPY(t, fullfile(other, 'spike_times.npy'));
            [pa, sa] = processed(ks); [pb, sb] = processed(other);
            tc.verifyNotEqual(sa.SortingID, sb.SortingID);
            [~, ia, ib] = intersect([pa.Units.TemplateID], [pb.Units.TemplateID]);
            tc.verifyEmpty(intersect(string({pa.Units(ia).UnitID}), string({pb.Units(ib).UnitID})));
        end

        function testUnitIDsIndependentOfMountPoint(tc)
            ks = make_mock_sorting(tc.Root, 5, 80, 20000, 60);
            s0 = makeSegment(ks, 'segment_0', 1);
            moved = fullfile(tc.Root, "mount2", "group", "sorter_output");
            mkdir(fileparts(moved)); copyfile(ks, moved);
            [pa, ~] = processed(s0);
            [pb, ~] = processed(fullfile(moved, 'segment_0'));
            tc.verifyEqual(string({pa.Units.UnitID}), string({pb.Units.UnitID}));
        end

        function testLegacySortingKey(tc)
            a = SpikeData.legacySortingKey("S:\lab\server01\data\Scope\user\Project\Plate1\well012\sorter_output\segment_3");
            b = SpikeData.legacySortingKey("/mnt/lab/data/Scope/user/Project/Plate1/well012/sorter_output/segment_0/");
            tc.verifyEqual(a, b);                                   % segments and mount points agree
            tc.verifyEqual(a, "Scope/user/Project/Plate1/well012/sorter_output");
            c = SpikeData.legacySortingKey("/x/Dataset/240101/CHIP01/Network/well000/sorter_output/qc_output");
            d = SpikeData.legacySortingKey("/x/Dataset/240109/CHIP01/Network/well000/sorter_output/qc_output");
            tc.verifyNotEqual(c, d);                                % other recording day
        end

    end
end

function seg = makeSegment(ks, name, half)
% copy of the sorting holding one half of its spikes, as a segment folder inside it
    seg = fullfile(ks, name);
    mkdir(seg);
    f = dir(ks); f = f(~[f.isdir]);
    for k = 1:numel(f), copyfile(fullfile(ks, f(k).name), seg); end
    t = readNPY(fullfile(seg, 'spike_times.npy')); u = readNPY(fullfile(seg, 'spike_templates.npy'));
    keep = (1:numel(t))' <= numel(t) / 2;
    if half == 2, keep = ~keep; end
    writeNPY(t(keep), fullfile(seg, 'spike_times.npy'));
    writeNPY(u(keep), fullfile(seg, 'spike_templates.npy'));
end

function [p, sd] = processed(path)
    md = struct('ChipID', "chip1", 'PlatingDate', "250101");
    [~, sd] = evalc('SpikeData.fromKilosort(path, md)');
    p = RecordingProcessor(sd, struct('QC', struct('N_Units', 1, 'RPV', 1)));
    evalc('p.runQC();');
end
