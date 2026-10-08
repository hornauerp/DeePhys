classdef test_detectNetworkBursts < matlab.unittest.TestCase
% TEST_DETECTNETWORKBURSTS  Joint ISI_N burst detection on synthetic networks with known bursts.
%
% Run: runtests('test_detectNetworkBursts')

    methods (Test)

        function testRecoversNetworkBursts(tc)
            [st, su, truth] = simulateNetworkBursts(Seed=1);
            [b, info] = detectNetworkBursts(st, su, 300);
            tc.verifyTrue(info.Gate);
            [prec, rec] = test_detectNetworkBursts.match(b, truth);
            tc.verifyGreaterThanOrEqual(prec, 0.9);
            tc.verifyGreaterThanOrEqual(rec, 0.9);
            % the regular, burst-blind units are the ones removed
            tc.verifyTrue(all(info.Units.Removed(truth.Tonic)));
            tc.verifyLessThan(mean(info.Units.Removed(~truth.Tonic)), 0.1);
        end

        function testOnsetsAndPeaks(tc)
            [st, su, truth] = simulateNetworkBursts(Seed=3);
            [b, info] = detectNetworkBursts(st, su, 300);
            err = nan(size(truth.T_start));
            for j = 1:numel(truth.T_start)
                h = find(b.T_start <= truth.T_end(j) & b.T_end >= truth.T_start(j), 1);
                if ~isempty(h), err(j) = abs(b.Onset(h) - truth.T_start(j)); end
            end
            tc.verifyLessThan(median(err, 'omitnan'), 0.02);
            tc.verifyTrue(all(b.T_peak >= b.Onset & b.T_peak <= b.Offset));
            % spans are first / last spikes: they contain the shape interval within +-0.1 s
            tc.verifyTrue(all(b.T_start <= b.Onset + 0.1 & b.T_end >= b.Offset - 0.1));
            % burst size counts the kept units' spikes between Onset and Offset
            kept = ismember(su, info.Units.Unit(~info.Units.Removed));
            in = kept & st >= b.Onset(1) & st <= b.Offset(1);
            tc.verifyEqual(b.NSpikes(1), sum(in));
            tc.verifyEqual(b.NUnits(1), numel(unique(su(in))));
        end

        function testNestedBurstsGiveTwoScales(tc)
            % three sub-events per network burst, < MinQuietGap apart (duty cycle 0.7: gaps <= 0.15 x
            % the burst duration, < 0.2 s even for the longest bursts): the bursts are the main
            % scale, the sub-events a second scale inside them
            [st, su, truth] = simulateNetworkBursts(SubBursts=3, BurstDuration=0.6, SubDutyCycle=0.7, Seed=1);
            [b, info] = detectNetworkBursts(st, su, 300);
            S = info.Scales;
            tc.verifyGreaterThanOrEqual(height(S), 2);
            tc.verifyEqual(numel(b.T_start), numel(truth.T_start));
            k = find(S.NBursts == numel(truth.Sub_start));
            tc.verifyNotEmpty(k);
            tc.verifyTrue(any(S.Inside(k) >= 0.8));
        end

        function testSmallEventsCountAsBursts(tc)
            % small synchronous events (15 of 60 units, two spikes each within 30 ms, far above
            % MinBurstUnits = 6) between large bursts: the operating point is the finest scale
            % that keeps the large bursts whole, so both kinds are found; coincidences of the
            % background are not (chance test)
            [st, su, truth] = simulateNetworkBursts(Seed=6);
            rs = RandStream('twister', 'Seed', 6);
            gap = [reshape(truth.T_end(1:end-1), 1, []); reshape(truth.T_start(2:end), 1, [])];
            gap = gap(:, diff(gap) > 2);
            ts = gap(1, 1:min(20, size(gap, 2))) + 1;               % 1 s after a large burst
            for e = ts
                u = randperm(rs, 54, 15)';                           % bursting (non-tonic) units
                st = [st; e + 0.03 * rand(rs, 30, 1)]; su = [su; repmat(u, 2, 1)]; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
            [b, info] = detectNetworkBursts(st, su, 300);
            tc.verifyTrue(info.Gate);
            hit = @(x) any(b.T_start <= x + 0.03 & b.T_end >= x);
            tc.verifyGreaterThanOrEqual(mean(arrayfun(hit, ts)), 0.8);              % small events
            tc.verifyGreaterThanOrEqual(mean(arrayfun(hit, (truth.T_start + truth.T_end) / 2)), 0.95);
            tc.verifyLessThanOrEqual(numel(b.T_start), numel(truth.T_start) + numel(ts) + 2);
            tc.verifyNotEmpty(info.Levels);
        end

        function testSuperburstsAreSplitIntoBursts(tc)
            % 4 s superbursts - all units at 8 Hz - holding five bursts (0.15 s at 60 Hz, every
            % 0.8 s), plus single bursts in a quiet background: the activity between the bursts of
            % a superburst is far below their peaks but far above the background, so ISI_N merges
            % each superburst; it is split into its bursts and reported as a superburst
            rs = RandStream('twister', 'Seed', 8); U = 60; T = 300;
            [st, su] = test_detectNetworkBursts.poissonUnits(rs, U, 0, T, 0.5);
            sbt = 10:20:290; single = 20:20:280; truth = [];
            for s0 = sbt
                [a, u] = test_detectNetworkBursts.poissonUnits(rs, U, s0, s0 + 4, 8); st = [st; a]; su = [su; u]; %#ok<AGROW>
                for k = 0:4
                    b0 = s0 + 0.2 + 0.8 * k; truth(end+1) = b0; %#ok<AGROW>
                    [a, u] = test_detectNetworkBursts.poissonUnits(rs, U, b0, b0 + 0.15, 60); st = [st; a]; su = [su; u]; %#ok<AGROW>
                end
            end
            for b0 = single
                truth(end+1) = b0; %#ok<AGROW>
                [a, u] = test_detectNetworkBursts.poissonUnits(rs, U, b0, b0 + 0.15, 60); st = [st; a]; su = [su; u]; %#ok<AGROW>
            end
            [st, o] = sort(st); su = su(o);
            [b, info] = detectNetworkBursts(st, su, T);
            found = arrayfun(@(x) any(b.T_start <= x + 0.15 & b.T_end >= x), truth);
            merged = arrayfun(@(k) sum(truth >= b.T_start(k) - 0.15 & truth <= b.T_end(k)) > 1, 1:numel(b.T_start));
            tc.verifyGreaterThanOrEqual(mean(found), 0.95);
            tc.verifyLessThanOrEqual(mean(merged), 0.05);
            tc.verifyEqual(size(info.Superbursts, 2), numel(sbt), 'AbsTol', 1);
            tc.verifyEqual(sum(b.Superburst > 0), 5 * numel(sbt), 'AbsTol', 4);
        end

        function testNoNetworkBurstsWithoutSynchrony(tc)
            % units burst on their own schedules: single-unit bursts, no network bursts
            for m = ["independent", "none"]
                [st, su] = simulateNetworkBursts(Mode=m, Seed=2);
                for tm = ["coupling", "none"]
                    [b, info] = detectNetworkBursts(st, su, 300, TonicMethod=tm);
                    tc.verifyFalse(info.Gate, m + " / " + tm);
                    tc.verifyEmpty(b.T_start, m + " / " + tm);
                end
            end
        end

        function testFixedParameterISIN(tc)
            % three spikes within 50 ms make a burst; a 0.48 s gap separates two bursts even
            % though no isolated spike lies between them
            t = [0 0.01 0.02 0.5 0.51 0.52 3]';
            b = detectNetworkBursts(t, ones(7, 1), 4, N=3, Threshold=0.05);
            tc.verifyEqual(b.T_start, [0 0.5], 'AbsTol', 1e-12);
            tc.verifyEqual(b.T_end, [0.02 0.52], 'AbsTol', 1e-12);
        end

        function testGeneratorIsReproducible(tc)
            [a, ua, ta] = simulateNetworkBursts(Seed=5);
            [b, ub, tb] = simulateNetworkBursts(Seed=5);
            tc.verifyEqual(a, b); tc.verifyEqual(ua, ub); tc.verifyEqual(ta, tb);
        end

    end

    methods (Static)
        function [t, u] = poissonUnits(rs, U, a, b, hz)
            % Poisson spikes of U units at hz in [a, b)
            t = []; u = [];
            for k = 1:U
                x = a + cumsum(-log(rand(rs, ceil(3 * hz * (b - a)) + 10, 1)) / hz);
                x = x(x < b); t = [t; x]; u = [u; k * ones(numel(x), 1)]; %#ok<AGROW>
            end
        end

        function [prec, rec] = match(b, truth)
            % a detected burst counts if it overlaps a true burst, and vice versa
            ov = @(s0, e0, s1, e1) any(s0 <= e1(:)' & e0 >= s1(:)', 2);
            prec = mean(ov(b.T_start(:), b.T_end(:), truth.T_start, truth.T_end));
            rec  = mean(ov(truth.T_start(:), truth.T_end(:), b.T_start, b.T_end));
        end
    end
end
