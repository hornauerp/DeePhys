function [feature_table, culture_array, values] = aggregateCultureFeatureTables(rg, level, metadata_field, values, tolerance, network_features, unit_features, normalization, useClustered)
            arguments
                rg RecordingGroup
                level string = "Recording"
                metadata_field string = "DIV"
                values {isnumeric} = [7,14,21,28] %refers to values of metadata_field; NaN to use the maximum number of available timepoints
                tolerance {isnumeric} = 1 %deviation from the actual age that will still be considered valid (e.g. age = 7 and tolerance = 1 considers DIVs 6-8)
                network_features string = "all"
                unit_features string = "all"
                normalization string = "baseline" % "baseline" (divided by first timepoint) or "scaled (scaled between 0 and 1)
                useClustered logical = false
            end

            rec_metadata = [rg.Recordings.Metadata];
            unique_values = unique([rec_metadata.(metadata_field)]);
            if isnan(values)
                values = unique_values;
            end
%             value_dev = abs(values - unique_values');
%             if tolerance == 0
%                 [r,c] = find(value_dev == 0);
%                 val_idx = unique(r);
%             else
%                 [r,c] = find(value_dev <= tolerance & value_dev > 0);
%                 val_idx = unique(min([r,c],[],2));
%             end

%             value_count = sum(value_dev<=tolerance);
%             if any(value_count == 0)
%                warning("DIV " + num2str(values(value_count == 0)) + " without corresponding cultures")
%             end
%             values = unique_values(val_idx);
            if isempty(values)
               error('No recordings found')
            end
            culture_table_array = {};
            culture_array = {};
            for iC = 1:length(rg.Cultures)
                rec_table_array = cell(1,length(values));
                culture_metadata = [rg.Cultures{iC}.Metadata];
                value = [culture_metadata.(metadata_field)];
                if length(value) >= length(values)
                    value_dev = abs(value - values');
                    value_count = sum(value_dev<=tolerance,1);
%                     value = value(value_count==1);
                    if sum(value_count==1) == length(values) % Need to find a way to handle two recordings falling within the tolerance window
                        sel_rec = rg.Cultures{iC}(value_count == 1);
                        sel_md = [sel_rec.Metadata];
                        [~, sort_idx] = sort([sel_md.(metadata_field)],'ascend');
                        sel_rec = sel_rec(sort_idx);
                        culture_array = [culture_array {sel_rec}];

                        for iR = 1:length(sel_rec)
                            if level == "Unit"
                                iR_table = MEArecording.getUnitFeatures(sel_rec(iR),unit_features);
                            elseif level == "Recording"
                                iR_table = getRecordingFeatures(sel_rec(iR), network_features, unit_features, useClustered);
                            else
                                error('Unknown level, select either "Unit" or "Recording"')
                            end
                            iR_table.Properties.VariableNames = iR_table.Properties.VariableNames + "_" + string(values(iR));
                            rec_table_array{iR} = iR_table;
                        end
                        if normalization == "baseline"
                            norm_mat = arrayfun(@(x) rec_table_array{x}.Variables./rec_table_array{1}.Variables,1:length(rec_table_array),'un',0);
                            norm_mat = [norm_mat{2:end}]; %2:end to omit initial 1
                            norm_mat(isnan(norm_mat)) = 0;
                            norm_mat(isinf(norm_mat)) = max(norm_mat(norm_mat < Inf),[],'all');

                            norm_table = [rec_table_array{2:end}]; %2:end to omit initial 1
                            norm_vars = norm_table.Properties.VariableNames;
                            culture_table = array2table(norm_mat,'VariableNames',norm_vars);
                            culture_table(:,var(culture_table.Variables) == 0) = [];
                            culture_table = [culture_table rec_table_array{1}(:,startsWith(string(rec_table_array{1}.Properties.VariableNames),"Waveform"))];
                            culture_table_array{iC} = culture_table;

                        elseif normalization == "scaled"
                            norm_mat = cellfun(@(x) x.Variables,rec_table_array,'un',0);
                            norm_mat = cat(3,norm_mat{:});
                            norm_mat = normalize(norm_mat,3,'range',[0 1]);
                            re_mat = reshape(norm_mat,size(norm_mat,1),[]);
                            norm_table = [rec_table_array{:}];
                            norm_vars = norm_table.Properties.VariableNames;
                            culture_table_array{iC} = array2table(re_mat,'VariableNames',norm_vars);
                        else
                            culture_table_array{iC} = [rec_table_array{:}];
                        end
                    else
                        continue
                    end
                else
                   continue
               end
            end
            clean_culture_tables = ~cellfun(@isempty, culture_table_array);
            culture_table_array = culture_table_array(clean_culture_tables);
            N_vars = cellfun(@width, culture_table_array);
            while length(unique(N_vars)) > 1
                [~,min_var_idx] = min(N_vars);
                min_vars = culture_table_array{min_var_idx}.Properties.VariableNames;
                culture_table_array = cellfun(@(x) x(:, matches(string(x.Properties.VariableNames), string(min_vars))),culture_table_array,'un',0);
                N_vars = cellfun(@width, culture_table_array);
                %warning('Not all features present in all cultures, reduced to shared features')
            end
            feature_table = vertcat(culture_table_array{:});
%             keep_idx = ~isnan(std(feature_table.Variables,'omitnan')); %Remove variables with 0 variance
%             feature_table = feature_table(:,keep_idx);
end

function [sparse_feature_mat, var_names, grouping_idx, separation_idx] = aggregateSparseFeatureTable(rg, network_features, unit_features, useClustered, grouping_var,...
                grouping_val, comparison_var, pooling_vals, tolerance)
            arguments
                rg RecordingGroup
                network_features string = "all"
                unit_features string = "all"
                useClustered logical = false
                grouping_var string = "DIV"
                grouping_val = nan
                comparison_var string= "Mutation"
                pooling_vals cell = {}
                tolerance = 1
            end
            rec_metadata = [rg.Recordings.Metadata];
            unique_values = unique([rec_metadata.(grouping_var)]);
            if isnan(grouping_val)
                grouping_val = unique_values;
            end
            test_table = getRecordingFeatures(rg.Recordings(r), network_features, unit_features, useClustered);
            num_fields = size(test_table,2);
            sparse_feature_mat = nan(length(rg.Cultures), length(grouping_val), num_fields);
            for i = 1:length(rg.Cultures)
                culture_metadata = [rg.Cultures{i}.Metadata];
                culture_vals = [culture_metadata.(grouping_var)];
                [vals, val_idx] = sort(culture_vals,'ascend');
                recordings = rg.Cultures{i}(val_idx);
                value_dev = abs(grouping_val - vals');
                value_count = find(sum(value_dev<=tolerance));
                for r = 1:length(recordings)
                    iR_table = getRecordingFeatures(recordings(r), network_features, unit_features, useClustered);
                    sparse_feature_mat(i,value_count(r),1:size(iR_table,2)) = iR_table.Variables;
                end
            end
            var_names = string(iR_table.Properties.VariableNames);
            [separation_idx, ~] = combineMetadataIndices(rg, rg.Cultures, comparison_var, pooling_vals); % e.g. mutation or treatment
            [grouping_idx, ~] = combineMetadataIndices(rg, rg.Cultures, grouping_var); %e.g. DIV or treatment concentrations

end
