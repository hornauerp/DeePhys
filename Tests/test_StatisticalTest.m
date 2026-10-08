classdef test_StatisticalTest < matlab.unittest.TestCase
% TEST_STATISTICALTEST  Group comparisons with units nested in cultures.
%
% Run: runtests('test_StatisticalTest')

    methods (Test)

        function testClustersRemovePseudoreplication(tc)
            % 4 + 4 cultures, 50 units each; cultures differ by a random offset but the
            % genotypes do not. Unit-level tests see 400 "independent" samples and are
            % often significant; culture medians are not.
            rng(41);
            n_sig_units = 0; n_sig_cult = 0;
            for rep = 1:20
                cult = repelem((1:8)', 50);
                geno = repelem(["WT"; "KO"], 200);
                v = repelem(randn(8, 1), 50) + 0.3 * randn(400, 1);
                r_u = StatisticalTest.compareGroups(v, geno, test="ranksum");
                r_c = StatisticalTest.compareGroups(v, geno, test="ranksum", Clusters=cult);
                n_sig_units = n_sig_units + r_u.significant;
                n_sig_cult  = n_sig_cult + r_c.significant;
                tc.verifyEqual(r_c.group_ns, [4 4]);
            end
            tc.verifyGreaterThan(n_sig_units, 8);    % far above the nominal 1 of 20
            tc.verifyLessThanOrEqual(n_sig_cult, 3);
        end

        function testClusterSpanningGroupsErrors(tc)
            v = randn(6, 1); g = ["a"; "a"; "b"; "b"; "a"; "b"]; c = [1; 1; 1; 2; 2; 2];
            tc.verifyError(@() StatisticalTest.compareGroups(v, g, Clusters=c), 'StatisticalTest:clusters');
        end

        function testAnovaWithClusters(tc)
            rng(42);
            cult = repelem((1:9)', 20);
            grp = repelem(["A"; "B"; "C"], 60);
            v = repelem([0; 0; 0; 0; 0; 0; 3; 3; 3], 20) + 0.2 * randn(180, 1);
            r = StatisticalTest.anova(v, grp, Clusters=cult, test="anova");   % 3 medians per group
            tc.verifyLessThan(r.p_omnibus, 0.05);
        end

    end
end
