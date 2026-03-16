function [mlm_result_array, p_G, p_GxT, features] = runMLM(rg, network_features, unit_features, feature_names, useClustered, grouping_var, comparison_var)
            arguments
                rg RecordingGroup
                network_features string = "all"
                unit_features string = "all"
                feature_names string = []
                useClustered logical = false
                grouping_var = "DIV"
                comparison_var = "Mutation"
            end

            [feat_mat, features, separation_val, grouping_val, subject_val] = prepareMLMinput(rg, network_features, unit_features, useClustered, grouping_var, comparison_var);

            if ~isempty(feature_names)
                feat_tbl = array2table(feat_mat,'VariableNames',features);
                sel_tbl = feat_tbl(:,feature_names);
                feat_mat = sel_tbl.Variables;
                features = feature_names;
            end

            mlm_result_array = cell(1,length(features));
            p_G = nan(size(features));
            p_GxT = nan(size(features));
            for f = 1:length(features)
                y = feat_mat(:,f); y(isinf(y)) = NaN;
                % Remove rows with NaN for this feature (fitlme cannot handle NaN)
                valid = ~isnan(y);
                tbl = table(y(valid), str2double(grouping_val(valid)),subject_val(valid),separation_val(valid),'VariableNames',["y","Week","Subject",comparison_var]);
                formula = ['y~ Week*', char(comparison_var), ' + (1|Subject)'];
                lme = fitlme(tbl,formula,'FitMethod','REML','DummyVarCoding','effects');
                mlm_result_array{f} = anova(lme,'DFMethod','satterthwaite');
                p_G(f) = mlm_result_array{f}.pValue(3) * length(features);
                p_GxT(f) = mlm_result_array{f}.pValue(4) * length(features);
            end

end

function [feat_mat, features, separation_val, grouping_val, subject_val] = prepareMLMinput(rg, network_features, unit_features, useClustered, grouping_var, comparison_var, pooling_vals)
            arguments
                rg RecordingGroup
                network_features string = "all"
                unit_features string = "all"
                useClustered logical = false
                grouping_var string = "DIV"
                comparison_var string = "Mutation"
                pooling_vals cell = {}
            end

            for r = 1:length(rg.Recordings)
                iR_table = getRecordingFeatures(rg.Recordings(r), network_features, unit_features, useClustered);
                feat_mat(r,:) = iR_table.Variables;
            end
            features = string(iR_table.Properties.VariableNames);
            [separation_idx, separation_labels] = combineMetadataIndices(rg, rg.Recordings, comparison_var, pooling_vals); % e.g. mutation or treatment
            [grouping_idx, grouping_labels] = combineMetadataIndices(rg, rg.Recordings, grouping_var); %e.g. DIV or treatment concentrations
            [subject_idx, subject_labels] = combineMetadataIndices(rg, rg.Recordings, "ChipID");
            separation_val = categorical(separation_labels(separation_idx));
            grouping_val = grouping_labels(grouping_idx);
            subject_val = categorical(subject_labels(subject_idx));
end
