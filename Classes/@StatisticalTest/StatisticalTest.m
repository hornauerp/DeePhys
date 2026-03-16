classdef StatisticalTest
% STATISTICALTEST  Framework for statistical comparisons on MEA data.
%
%   Provides paired/unpaired group comparisons, mixed-effects models for
%   longitudinal data, and multiple comparison correction. All tests return
%   standardized result structs.
%
%   USAGE:
%     % Compare feature across two groups:
%     result = StatisticalTest.compareGroups(feat_values, group_labels);
%
%     % Compare feature across multiple groups with correction:
%     result = StatisticalTest.anova(feat_values, group_labels, 'correction', 'bonferroni');
%
%     % Paired comparison (pre vs post treatment within cultures):
%     result = StatisticalTest.pairedComparison(pre_values, post_values, culture_ids);
%
%     % Mixed-effects model:
%     result = StatisticalTest.mixedEffects(feat_table, 'Feature ~ DIV * Mutation + (1|CultureID)');

    methods (Static)

        function result = compareGroups(values, group_labels, options)
        % COMPAREGROUPS  Two-group comparison with automatic test selection.
        %
        %   Uses Wilcoxon rank-sum (Mann-Whitney U) by default. Optionally
        %   uses t-test if data passes normality check.
        %
        % INPUTS:
        %   values       - (N x 1) numeric values
        %   group_labels - (N x 1) group assignments (exactly 2 unique values)
        %   options      - name-value: 'test' ("auto"|"ranksum"|"ttest"), 'alpha' (0.05)
        %
        % OUTPUT:
        %   result - struct with: p, stat, test_name, effect_size (Cohen's d),
        %            group_means, group_medians, group_ns, significant
            arguments
                values (:,1) double
                group_labels
                options.test string = "auto"
                options.alpha (1,1) double = 0.05
            end

            groups = unique(group_labels);
            assert(length(groups) == 2, 'StatisticalTest:compareGroups requires exactly 2 groups. Use anova for >2.');

            g1 = values(group_labels == groups(1));
            g2 = values(group_labels == groups(2));

            % Automatic test selection
            if options.test == "auto"
                if length(g1) >= 8 && length(g2) >= 8
                    [~, p_norm1] = kstest((g1 - mean(g1)) / std(g1));
                    [~, p_norm2] = kstest((g2 - mean(g2)) / std(g2));
                    if p_norm1 > 0.05 && p_norm2 > 0.05
                        test_name = "ttest";
                    else
                        test_name = "ranksum";
                    end
                else
                    test_name = "ranksum";
                end
            else
                test_name = options.test;
            end

            switch test_name
                case "ttest"
                    [~, p, ~, stats] = ttest2(g1, g2);
                    stat = stats.tstat;
                case "ranksum"
                    [p, ~, stats] = ranksum(g1, g2);
                    if isfield(stats, 'zval')
                        stat = stats.zval;
                    else
                        stat = stats.ranksum;
                    end
            end

            % Effect size (Cohen's d)
            pooled_sd = sqrt(((length(g1)-1)*var(g1) + (length(g2)-1)*var(g2)) / (length(g1)+length(g2)-2));
            if pooled_sd > 0
                cohens_d = (mean(g1) - mean(g2)) / pooled_sd;
            else
                cohens_d = 0;
            end

            result.p = p;
            result.stat = stat;
            result.test_name = test_name;
            result.effect_size = cohens_d;
            result.group_means = [mean(g1), mean(g2)];
            result.group_medians = [median(g1), median(g2)];
            result.group_ns = [length(g1), length(g2)];
            result.group_labels = groups;
            result.significant = p < options.alpha;
            result.alpha = options.alpha;
        end

        function result = anova(values, group_labels, options)
        % ANOVA  One-way ANOVA or Kruskal-Wallis with post-hoc correction.
        %
        % INPUTS:
        %   values       - (N x 1) numeric values
        %   group_labels - (N x 1) group assignments (>=2 groups)
        %   options      - name-value: 'correction' ("bonferroni"|"holm"|"fdr"),
        %                  'alpha' (0.05), 'test' ("auto"|"anova"|"kruskalwallis")
        %
        % OUTPUT:
        %   result - struct with: p_omnibus, stat, test_name, posthoc (table),
        %            group_means, significant
            arguments
                values (:,1) double
                group_labels
                options.correction string = "bonferroni"
                options.alpha (1,1) double = 0.05
                options.test string = "auto"
            end

            groups = unique(group_labels);
            n_groups = length(groups);

            % Test selection
            if options.test == "auto"
                % Check normality of each group
                all_normal = true;
                for g = 1:n_groups
                    gv = values(group_labels == groups(g));
                    if length(gv) >= 8
                        [~, p_ks] = kstest((gv - mean(gv)) / std(gv));
                        if p_ks <= 0.05
                            all_normal = false;
                            break;
                        end
                    else
                        all_normal = false;
                        break;
                    end
                end
                test_name = string(all_normal * "anova" + ~all_normal * "kruskalwallis");
                if test_name == "0"
                    test_name = "kruskalwallis";
                elseif test_name == "1"
                    test_name = "anova";
                end
            else
                test_name = options.test;
            end

            switch test_name
                case "anova"
                    [p_omnibus, tbl, stats] = anova1(values, group_labels, 'off');
                    stat = tbl{2, 5}; % F statistic
                case "kruskalwallis"
                    [p_omnibus, tbl, stats] = kruskalwallis(values, group_labels, 'off');
                    stat = tbl{2, 5}; % chi-sq statistic
            end

            % Post-hoc pairwise comparisons
            n_pairs = n_groups * (n_groups - 1) / 2;
            posthoc = table('Size', [n_pairs, 5], ...
                'VariableTypes', {'string', 'string', 'double', 'double', 'logical'}, ...
                'VariableNames', {'Group1', 'Group2', 'p_raw', 'p_corrected', 'significant'});
            idx = 0;
            raw_ps = zeros(n_pairs, 1);
            for i = 1:n_groups
                for j = (i+1):n_groups
                    idx = idx + 1;
                    g1 = values(group_labels == groups(i));
                    g2 = values(group_labels == groups(j));
                    raw_ps(idx) = ranksum(g1, g2);
                    posthoc.Group1(idx) = string(groups(i));
                    posthoc.Group2(idx) = string(groups(j));
                    posthoc.p_raw(idx) = raw_ps(idx);
                end
            end

            % Multiple comparison correction
            switch options.correction
                case "bonferroni"
                    corrected_ps = min(raw_ps * n_pairs, 1);
                case "holm"
                    corrected_ps = StatisticalTest.holmCorrection(raw_ps);
                case "fdr"
                    corrected_ps = StatisticalTest.fdrCorrection(raw_ps);
                otherwise
                    corrected_ps = raw_ps;
            end
            posthoc.p_corrected = corrected_ps;
            posthoc.significant = corrected_ps < options.alpha;

            % Group means
            group_means = arrayfun(@(g) mean(values(group_labels == g)), groups);

            result.p_omnibus = p_omnibus;
            result.stat = stat;
            result.test_name = test_name;
            result.posthoc = posthoc;
            result.correction = options.correction;
            result.group_means = group_means;
            result.group_labels = groups;
            result.significant = p_omnibus < options.alpha;
            result.alpha = options.alpha;
        end

        function result = pairedComparison(pre_values, post_values, pair_ids, options)
        % PAIREDCOMPARISON  Paired test for pre/post within the same group (e.g., culture).
        %
        % INPUTS:
        %   pre_values  - (N x 1) values before treatment
        %   post_values - (N x 1) values after treatment
        %   pair_ids    - (N x 1) pairing variable (e.g., culture ID)
        %   options     - name-value: 'test' ("auto"|"signrank"|"ttest"), 'alpha' (0.05)
            arguments
                pre_values (:,1) double
                post_values (:,1) double
                pair_ids = []
                options.test string = "auto"
                options.alpha (1,1) double = 0.05
            end

            assert(length(pre_values) == length(post_values), 'Pre and post must have same length.');
            diffs = post_values - pre_values;

            if options.test == "auto"
                if length(diffs) >= 8
                    [~, p_norm] = kstest((diffs - mean(diffs)) / std(diffs));
                    test_name = string(p_norm > 0.05) * "ttest" + string(p_norm <= 0.05) * "signrank";
                    if test_name == "1" || test_name == "0"
                        test_name = "signrank"; % fallback
                    end
                else
                    test_name = "signrank";
                end
            else
                test_name = options.test;
            end

            switch test_name
                case "ttest"
                    [~, p, ~, stats] = ttest(pre_values, post_values);
                    stat = stats.tstat;
                case "signrank"
                    [p, ~, stats] = signrank(pre_values, post_values);
                    if isfield(stats, 'zval')
                        stat = stats.zval;
                    else
                        stat = stats.signedrank;
                    end
            end

            % Effect size (paired Cohen's d)
            if std(diffs) > 0
                cohens_d = mean(diffs) / std(diffs);
            else
                cohens_d = 0;
            end

            result.p = p;
            result.stat = stat;
            result.test_name = test_name;
            result.effect_size = cohens_d;
            result.mean_diff = mean(diffs);
            result.median_diff = median(diffs);
            result.n_pairs = length(diffs);
            result.significant = p < options.alpha;
            result.alpha = options.alpha;
        end

        function result = mixedEffects(data_table, formula)
        % MIXEDEFFECTS  Fit a linear mixed-effects model.
        %
        %   Wrapper around MATLAB's fitlme for longitudinal/hierarchical data.
        %
        % INPUTS:
        %   data_table - table with columns matching the formula variables
        %   formula    - Wilkinson notation string, e.g. 'Feature ~ DIV * Mutation + (1|CultureID)'
        %
        % OUTPUT:
        %   result - struct with: lme (fitted model), anova_table, coefficients, R2
            arguments
                data_table table
                formula string
            end

            lme = fitlme(data_table, formula);
            result.lme = lme;
            result.anova_table = anova(lme);
            result.coefficients = lme.Coefficients;
            result.R2_conditional = lme.Rsquared.Ordinary;
            result.R2_marginal = lme.Rsquared.Adjusted;
            result.formula = formula;
        end

        function corrected = holmCorrection(raw_ps)
        % HOLMCORRECTION  Holm-Bonferroni step-down correction.
            n = length(raw_ps);
            [sorted_p, sort_idx] = sort(raw_ps);
            corrected = zeros(n, 1);
            for i = 1:n
                corrected(sort_idx(i)) = min(sorted_p(i) * (n - i + 1), 1);
            end
            % Enforce monotonicity
            [~, sort_idx2] = sort(raw_ps);
            for i = 2:n
                corrected(sort_idx2(i)) = max(corrected(sort_idx2(i)), corrected(sort_idx2(i-1)));
            end
        end

        function corrected = fdrCorrection(raw_ps)
        % FDRCORRECTION  Benjamini-Hochberg false discovery rate correction.
            n = length(raw_ps);
            [sorted_p, sort_idx] = sort(raw_ps);
            corrected = zeros(n, 1);
            for i = 1:n
                corrected(sort_idx(i)) = min(sorted_p(i) * n / i, 1);
            end
            % Enforce monotonicity (from largest to smallest)
            [~, rev_sort] = sort(raw_ps, 'descend');
            for i = 2:n
                corrected(rev_sort(i)) = min(corrected(rev_sort(i)), corrected(rev_sort(i-1)));
            end
        end

    end
end
