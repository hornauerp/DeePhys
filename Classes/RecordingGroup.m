classdef RecordingGroup < handle
    
    %%%%%%%TODO%%%%%%%
    % 
    %%%%%%%%%%%%%%%%%%
    
    properties
        Recordings
        Units
        Cultures
        Parameters %Inclusion // exclusion criteria
        DimensionalityReduction
        Clustering
        Classification
        AppliedClassification
        Regression
    end
    
    methods (Hidden)
        function rg = RecordingGroup(recording_array, parameters)
            arguments
                recording_array MEArecording
                parameters struct = struct();
            end
            
            rg.parseParameters(parameters);
            rg.harmonizeMetadata(recording_array);
            keep_recordings = rg.filterRecordingArray([recording_array.Metadata]);
            rg.Recordings = recording_array(keep_recordings);
            rg.Units = [rg.Recordings.Units];
            if isfield([rg.Recordings.Metadata],'ChipID')
                rg.Cultures = rg.groupCultures();
            end
            fprintf('Initialized RecordingGroup with %i recordings and %i units\n',...
                length(rg.Recordings),length(rg.Units))
        end
        
        function parseParameters(rg, parameters)
            rg.Parameters = rg.returnDefaultParams();
            parameter_fields = fieldnames(parameters);
            if isempty(parameter_fields)
                warning("No parameters provided, using default ones")
            else
                for pf = 1:length(parameter_fields)
                    if ismember(parameter_fields{pf},fieldnames(rg.Parameters))
                        subfields = fieldnames(parameters.(parameter_fields{pf}));
                        for sf = 1:length(subfields)
                            if ismember(subfields{sf},fieldnames(rg.Parameters.(parameter_fields{pf})))
                                rg.Parameters.(parameter_fields{pf}).(subfields{sf}) = parameters.(parameter_fields{pf}).(subfields{sf});
                            else
                                warning("Unrecognized parameter field: %s.%s",parameter_fields{pf}{:},subfields{sf}{:})
                            end
                        end
                    else
                        warning("Unrecognized parameter field: %s",pf{:})
                    end
                end
                disp("Imported custom parameters")
            end
        end
        
        function filtered_idx = filterRecordingArray(rg, metadata_array, inclusion, exclusion)
            arguments
                rg RecordingGroup
                metadata_array
                inclusion = rg.Parameters.Selection.Inclusion;
                exclusion = rg.Parameters.Selection.Exclusion;
            end
            
            % Input as cell arrays
            for i = 1:length(inclusion)
                if isnumeric(inclusion{i}{2}) || isstring((inclusion{i}{2}))
                    value_vector = [metadata_array.(inclusion{i}{1})];
                    idx = ismember(value_vector,[inclusion{i}{2:end}]);
                elseif iscell(inclusion{i}{2})
                    value_vector = {metadata_array.(inclusion{i}{1})};
                    idx = ismember(value_vector,inclusion{i}{2:end});
                else
                    value_vector = [metadata_array.(inclusion{i}{1})];
                    idx = ismember(value_vector,inclusion{i}(2:end));
                end
                if ~exist('inclusion_idx','var')
                    inclusion_idx = idx;
                else
                    inclusion_idx = idx & inclusion_idx;
                end
            end
            if ~exist('inclusion_idx','var')
               inclusion_idx = ones(1,length(metadata_array)); 
            end
            
            for i = 1:length(exclusion)
                if isnumeric(exclusion{i}{2}) || isstring((exclusion{i}{2}))
                    value_vector = [metadata_array.(exclusion{i}{1})];
                    idx = ismember(value_vector,[exclusion{i}{2:end}]);
                else
                    value_vector = {metadata_array.(exclusion{i}{1})};
                    idx = ismember(value_vector,exclusion{i}(2:end));
                end
                if ~exist('exclusion_idx','var')
                    exclusion_idx = idx;
                else
                    exclusion_idx = idx | exclusion_idx;
                end
            end
            if ~exist('exclusion_idx','var')
               exclusion_idx = zeros(1,length(metadata_array)); 
            end
            
            filtered_idx = inclusion_idx & ~exclusion_idx;
        end
        
        function harmonizeMetadata(rg, recording_array)
            md_fields = arrayfun(@(x) fields(x.Metadata),recording_array,'un',0);
            unique_fields = unique(string(vertcat(md_fields{:})));
            for r = 1:length(recording_array)
                missing_field = unique_fields(~matches(unique_fields,md_fields{r}));
                for mf = 1:length(missing_field)
                    recording_array(r).Metadata.(missing_field(mf)) = [];
                end
            end
        end
        
        function cultures = groupCultures(rg, grouping_var)
            arguments
               rg RecordingGroup
               grouping_var = [] %Needed if recordings from the same culture need to be grouped at different time points (e.g. two different dose responses)
            end
            cultures = [];
            metadata = [rg.Recordings.Metadata];
            if ~isfield(metadata,"ChipID")
                warning('No chip IDs provided, cannot group recordings')
                return
            end
            chip_ids = unique([metadata.ChipID]);
            if isempty(chip_ids)
                warning('No chip IDs provided, cannot group recordings')
                return
            end
            for id = chip_ids
                inclusion = {{'ChipID',id}};
                chip_idx = rg.filterRecordingArray(metadata, inclusion);
                chip_array = rg.Recordings(chip_idx);
                chip_metadata = [chip_array.Metadata];
                if ~isfield(metadata,"PlatingDate")
                    warning('No plating date provided, cannot differentiate batches')
                    return
                end
                chip_plating_dates = unique([chip_metadata.PlatingDate]);
                for pd = chip_plating_dates
                   inclusion = {{'PlatingDate',pd}};
                   plating_idx = rg.filterRecordingArray(chip_metadata, inclusion);
                   culture_array = chip_array(plating_idx);
                   if ~isempty(grouping_var)
                       group_metadata = [culture_array.Metadata];
                       group_values = unique([group_metadata.(grouping_var)]);
                       for gv = group_values
                           inclusion = {{grouping_var, gv}};
                           group_idx = rg.filterRecordingArray(group_metadata, inclusion);
                           group_array = culture_array(group_idx);
                           cultures = [cultures {group_array}];
                       end
                   else
                       cultures = [cultures {culture_array}];
                   end
                end
            end
        end
    end
    
    methods (Static)
       function defaultParams = returnDefaultParams()
           %Input to the filterRecordingArray function
           defaultParams.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
           defaultParams.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
       end
       
       function aligned_wf = alignWaveforms(unit_array, interp_factor)
           arguments
               unit_array Unit
               interp_factor double = 10
           end
           ref_wf = [unit_array.ReferenceWaveform];
           ref_wf = ref_wf(sum(ref_wf,2)~=0,:);
           [~,i] = min(ref_wf,[],1);
           peak_idx = mean(i);
           max_offset = round(peak_idx/2);
           x = max_offset:size(ref_wf,1)+max_offset-1;
           % xq = (1:size(ref_wf,1)+2*max_offset) * interp_factor;
           buffer_wf = size(ref_wf,1)+2*max_offset;
           xq = linspace(1,buffer_wf,buffer_wf*interp_factor);

           % interp_wf = interp1(x,ref_wf,xq,"linear",'extrap');
           interp_wf = interp1(x,ref_wf,xq,"makima");
           interp_wf = interp_wf((max_offset*interp_factor)+1:((buffer_wf-max_offset-1)*interp_factor),:);
           interp_wf = interp_wf./max(abs(interp_wf));
           rm_idx = find(ceil(abs(peak_idx - i)) >= max_offset);
           [~, minIndices] = min(interp_wf);
           shifts = round(peak_idx*interp_factor) - minIndices;

           aligned_wf = zeros(size(interp_wf));
           for j = 1:size(interp_wf, 2)
               if j == rm_idx
                   aligned_wf(:, j) = interp_wf(:, j);
               else
                   aligned_wf(:, j) = circshift(interp_wf(:, j), [shifts(j), 0]);
               end
           end
           shifts(rm_idx) = [];
           % aligned_wf = aligned_wf(max(shifts):end+min(shifts),:);
           aligned_wf = aligned_wf((5*interp_factor):(end-5*interp_factor),:);
       end
       
       function [clf,train_acc] = create_classifier(X_train,Y_train,alg,N_hyper)
           if N_hyper > 0
               switch alg
                   case 'svm'
                       clf = fitcsvm(X_train,Y_train,'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
                           struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                           'Verbose',0));
                   case 'cnb'
                       clf = fitcnb(X_train,Y_train,'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
                           struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                           'Verbose',0));
                   case 'knn'
                       clf = fitcknn(X_train,Y_train,'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',...
                           struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                           'Verbose',0));
                   case 'rf'
                       hyperparams = {'NumLearningCycles','MinLeafSize','MaxNumSplits','SplitCriterion','NumVariablesToSample'};
                       t = templateTree('Reproducible',true);
                       clf = fitcensemble(X_train,Y_train,'Method','Bag','OptimizeHyperparameters',hyperparams,'Learners',t, ...
                           'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',N_hyper,'ShowPlots',false,...
                           'Verbose',0));
               end
               train_acc = clf.HyperparameterOptimizationResults.MinObjective;
           else
               switch alg
                   case 'svm'
                       clf = fitcsvm(X_train,Y_train);
                   case 'cnb'
                       clf = fitcnb(X_train,Y_train);
                   case 'knn'
                       clf = fitcknn(X_train,Y_train);
                   case 'rf'
                       t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all','Reproducible',true);
                       clf = fitcensemble(X_train,Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t,'Options',statset("UseParallel",true));
               end
               train_acc = 1-resubLoss(clf,'LossFun','classiferror');
           end
       end
       
       function [Y_train, Y_test, train_idx, test_idx] = cv_split(Y, cv, k)
          arguments
             Y
             cv cvpartition
             k double
          end
          
           if iscell(Y)
               Y_train = [Y{cv.training(k)}];
               Y_test = [Y{cv.test(k)}];
               k_train = cv.training(k);
               train_idx = arrayfun(@(x) ones(size(Y{x})) * k_train(x),1:length(Y),'un',0);
               train_idx = [train_idx{:}];
               test_idx = ~train_idx;
               
           else
               Y_train = Y(cv.training(k));
               Y_test = Y(cv.test(k));
               train_idx = cv.training(k);
               test_idx = cv.test(k);
           end
       end
       
       function [new_group_idx,new_group_labels] = poolMetadataValues(group_idx, group_labels, classification_val)
           arguments
               group_idx double %Original group idx, corresponding to group_labels
               group_labels string
               classification_val string %Value(s) that are to be pooled (against)
           end
           clf_group_idx = find(contains(group_labels, classification_val));
           new_group_idx = ismember(group_idx,clf_group_idx) * 1;%(group_idx == clf_group_idx) * 1;
           new_group_labels(1) = join(group_labels(clf_group_idx),'/');
           clf_group_idx = find(~ismember(group_labels, classification_val));
           new_group_idx(new_group_idx == 0) = 2;
           new_group_labels(2) = join(group_labels(clf_group_idx),'/');
       end
       
       function [reduction, cmap] = plot_cluster_outlines(reduction, cluster_idx, ax, plot_centroid, nodeSz, mapSz, sigma, cmap)
           arguments
               reduction {isnumeric}
               cluster_idx (1,:) {isnumeric}
               ax = axes;
               plot_centroid = false
               nodeSz (1,1) {isnumeric} = 10
               mapSz {isnumeric} = 500
               sigma (1,1) {isnumeric} = mapSz/6
               cmap (:,3) {isnumeric} = othercolor('Spectral6',max(cluster_idx)); %othercolor("Set29",4);%
           end
           
           buffer_ratio = 0.5;
           reduction = reduction - min(min(reduction)) * (1 + buffer_ratio);
           
           xrng = ...                      % range of x coordinates (with buffer)
               [min(reduction(:,1)) - buffer_ratio*range(reduction(:,1)),max(reduction(:,1)) + buffer_ratio*range(reduction(:,1))];
           yrng = ...                      % range of y coordinates (with buffer)
               [min(reduction(:,2)) - buffer_ratio*range(reduction(:,2)),max(reduction(:,2)) + buffer_ratio*range(reduction(:,2))];
           
           grid = zeros(mapSz,mapSz,max(cluster_idx));
           for i = 1:max(cluster_idx)                  % loop over all categories
               gridtemp = zeros(mapSz);        % temporary grid
               idx = cluster_idx == i;                 % get current category
               coori = reduction(idx,:);            % coordinates
               x = interp1(xrng,[1,mapSz],coori(:,1)); % map x coordinates to cells
               y = interp1(yrng,[1,mapSz],coori(:,2)); % same for y
               rx = round(x);                  % round to whole #
               ry = round(y);                  % ditto
               pts = (rx - 1)*mapSz + ry;
               gridtemp(pts) = 1;
               grid(:,:,i) = gridtemp;
           end
           reduction(:,1) = interp1(xrng,[1,mapSz],reduction(:,1));
           reduction(:,2) = interp1(yrng,[1,mapSz],reduction(:,2));
           
           sz = ceil(sqrt(mapSz));
           x = linspace(-sz,sz);
           y = linspace(-sz,sz);
           [mx,my] = meshgrid(x,y);
           kernel = exp(-(mx.^2 + my.^2)./(2*sigma.^2));
           kernel = kernel/sum(kernel(:));
           
           tot = zeros(mapSz,mapSz,3);
           colall = zeros(mapSz,mapSz,3,max(cluster_idx));
           for i = 1:max(cluster_idx)
               gridsmooth = conv2(grid(:,:,i),kernel,'same');
               prob = gridsmooth/max(gridsmooth(:));
               tot(:,:,i) = prob;
               
               P = bsxfun(@minus,ones(mapSz,mapSz,3),prob);
               C = ones(mapSz,mapSz,3);
               for k = 1:3
                   C(:,:,k) = cmap(i,k)*prob;
               end
               col = P + C;
               colall(:,:,:,i) = col;
           end
           
           [~,idx] = max(tot,[],3);
           cc = ones(mapSz,mapSz,3);
           for i = 1:max(cluster_idx)
               mask = idx == i;
               [u,v] = find(mask);
               for j = 1:length(u)
                   cc(u(j),v(j),:) = colall(u(j),v(j),:,i);
               end
           end
           
           for i = 1:3
               cc(:,:,i) = conv2(cc(:,:,i),kernel,'same');
           end
           %            prct = round(mapSz*0.1);
           %            cc(1:prct,:,:) = 1;
           %            cc(:,1:prct,:) = 1;
           %            cc(:,(end - prct + 1):end,:) = 1;
           %            cc((end - prct + 1):end,:,:) = 1;
           
           hold(ax,'on');
           imagesc(cc)
           for i = 1:max(cluster_idx)
               idx = cluster_idx == i;
               coori = reduction(idx,:);
               if plot_centroid
                   s = scatter(coori(:,1),coori(:,2),nodeSz,repmat(cmap(i,:),sum(idx),1),'o','filled', 'MarkerFaceAlpha',0.8, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.25);
                   s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                   scatter(mean(coori(:,1)),mean(coori(:,2)),nodeSz*2,cmap(i,:),'o','filled', 'MarkerEdgeColor','k')
               else
                   scatter(coori(:,1),coori(:,2),nodeSz,repmat(cmap(i,:),sum(idx),1),'o','filled', 'MarkerEdgeColor','k');
               end
           end
           xlim([mapSz*0.1 mapSz*0.9])
           ylim([mapSz*0.1 mapSz*0.9])
           axis off
           %            axis image square;
       end
       
    end
    
    methods 
        
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
                y = feat_mat(:,f); y(isnan(y) | isinf(y)) = 0;
                tbl = table(y, str2double(grouping_val),subject_val,separation_val,'VariableNames',["y","Week","Subject",comparison_var]);
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
        
        function runANOVA(rg, level, metadata_field, values, tolerance, network_features, unit_features, normalization, useClustered)
           arguments
               rg RecordingGroup
               level string = "Recording"
               metadata_field string = "DIV"
               values {isnumeric} = [7,14,21,28] %refers to values of metadata_field // NaN to use the maximum number of available timepoints
               tolerance {isnumeric} = 1 %deviation from the actual age that will still be considered valid (e.g. age = 7 and tolerance = 1 considers DIVs 6-8)
               network_features string = "all"
               unit_features string = "all"
               normalization string = "baseline" % "baseline" (divided by first timepoint) or "scaled (scaled between 0 and 1)
               useClustered logical = false
           end
           
           [feature_table, culture_array, values] = aggregateCultureFeatureTables(rg, level, metadata_field, values, tolerance, network_features, unit_features, normalization, useClustered);
           %%% TO DO %%%
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
        
        function [feature_table, culture_array, values] = aggregateCultureFeatureTables(rg, level, metadata_field, values, tolerance, network_features, unit_features, normalization, useClustered)
            arguments
                rg RecordingGroup
                level string = "Recording"
                metadata_field string = "DIV"
                values {isnumeric} = [7,14,21,28] %refers to values of metadata_field // NaN to use the maximum number of available timepoints
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
        
        function [norm_train_data, norm_test_data] = normalizeByGroup(rg, feat_mat, object_group, grouping_var, train_idx, test_idx)
            arguments
                rg RecordingGroup
                feat_mat {isnumeric}
                object_group %Units/recordings/cultures corresponding to the rows of feat_mat
                grouping_var string = "PlatingDate" %Has to correspond to a MEArecording metadata field
                train_idx (1,:) logical = logical(ones(1,size(feat_mat,1))) %Default to unsupervised // should be binary and length of object_group
                test_idx (1,:) logical = logical(zeros(1,size(feat_mat,1)))
            end
            
            train_idx = logical(train_idx);
            test_idx = logical(test_idx);
            
            % switch class(object_group)
            %     case 'Unit'
            %         recordings = [object_group.MEArecording];
            %         metadata = [recordings.Metadata];
            %         metadata = [metadata.(grouping_var)];
            % 
            %     case 'MEArecording'
            %         metadata = [object_group.Metadata];
            %         metadata = [metadata.(grouping_var)];
            % 
            %     case 'cell'
            %         metadata = cellfun(@(x) string(x(1).Metadata.(grouping_var)),object_group);
            % end
            [iG, G] = rg.combineMetadataIndices(object_group, grouping_var);
            g_idx = unique(iG);
            % [iG,G] = findgroups(metadata);
            for g = 1:length(g_idx)
                iBatch_train = (iG == g_idx(g));% & train_idx';
                [norm_mat, batch_mean, batch_sd] = normalize(feat_mat(iBatch_train,:));
                feat_mat(iBatch_train,:) = norm_mat;%./max(abs(norm_mat)); %Check
                % if any(test_idx)
                %     iBatch_test = iG == g & test_idx';
                %     feat_mat(iBatch_test,:) = normalize(feat_mat(iBatch_test,:),'center',batch_mean,'scale',batch_sd);
                % end
            end
            if length(G) > 1 %Only normalize again if more than 1 group exists
                [norm_train_data, train_mean, train_sd] = normalize(feat_mat(train_idx,:));
                norm_test_data = normalize(feat_mat(test_idx,:),'center',train_mean,'scale',train_sd);
            else
                norm_train_data = feat_mat(train_idx,:);
                norm_test_data = feat_mat(test_idx,:);
            end
            % scale_factor = max(abs(norm_train_data));
            % norm_train_data = norm_train_data./scale_factor;
            % norm_test_data = norm_test_data./scale_factor;
        end
        
        function [norm_train, norm_test, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx)
            arguments
                rg RecordingGroup
                input_table table
                object_group %Units/recordings/cultures corresponding to the rows of feat_mat
                normalization_var string = "PlatingDate" %Has to correspond to a MEArecording metadata field
                train_idx (1,:) logical = logical(ones(1,size(input_table,1))) %Default to unsupervised // should be logical and length of object_group
                test_idx (1,:) logical = logical(zeros(1,size(input_table,1)))
            end
            
            input_mat = input_table.Variables;
            feature_names = string(input_table.Properties.VariableNames);
            input_mat(isnan(input_mat)) = 0; %Handle rare case where NaN appears
            if isempty(normalization_var)
                [norm_train, mean_train, sd_train] = normalize(input_mat(train_idx,:));
                nan_idx = any(isnan(norm_train));
                norm_train(:,nan_idx) = [];
                feature_names(nan_idx) = [];
                scale_factor = max(abs(norm_train)); %Make sure we go below the scale of UMAP
                norm_train = norm_train./scale_factor; %Scale data
                if sum(test_idx) > 0
                    
                    norm_test = normalize(input_mat(test_idx,:),'center',mean_train,'scale',sd_train);
                    norm_test(:,nan_idx) = [];
                    norm_test = norm_test./scale_factor;
                else
                    norm_test = [];
                end
            else
                [norm_train, norm_test] = normalizeByGroup(rg, input_mat, object_group, normalization_var, train_idx, test_idx); %Normalize
                nan_idx = any(isnan(norm_train)) | any(isnan(norm_test),1);
                norm_train(:,nan_idx) = [];
                feature_names(nan_idx) = [];
                norm_train = norm_train./max(abs(norm_train)); %Scale data
                if sum(test_idx) > 0
                    norm_test(:,nan_idx) = [];%Remove peak NaNs
                    norm_test = norm_test./max(abs(norm_train)); %Scale data
                else
                    norm_test = [];
                end
            end
            
        end
        
        function [reduction, input_table] = reduceDimensionality(rg, level, method, n_dims, unit_features, network_features, feature_names, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance)
            arguments
               rg RecordingGroup
               level string = "Unit" %Unit, Recording or Culture level
               method string = "UMAP" %UMAP, PCA
               n_dims (1,1) {isnumeric} = 2 %Number of output dimensions
               unit_features string = ["ReferenceWaveform","ActivityFeatures"]
               network_features string = "all"
               feature_names string = []
               useClustered logical = false
               normalization_var string = "PlatingDate"
               grouping_var string = [] %Has to correspond to a MEArecording metadata field
               grouping_values {isnumeric} = [7,14,21,28] %Only relevant for DimensionalityReduction on the culture level
               normalization string = [] %"scaled" or "baseline" or []
               tolerance {isnumeric} = 1 %Gives tolerance for culture selection by age (e.g. age=7 tolerance=1 allows DIVs 6-8)
            end
            
%             fprintf('Performing %s dimensionality reduction on %ss and normalizing by %s\n', method, level, grouping_var)
            if isempty(grouping_var) %Get group indices for stratification
                object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            end
            
            if isempty(grouping_var)
                if level == "Unit"
                    input_table = MEArecording.getUnitFeatures(object_group,unit_features);
                    object_group = [rg.Units];
                elseif level == "Recording"
                    input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                else 
                    error('Unknown level')
                end
            else
                if level == "Unit"
                    object_group = cellfun(@(x) [x(1).Units],object_group,'un',0);
                    object_group = [object_group{:}];
                end
            end
            
            if ~isempty(feature_names)
                vars = string(input_table.Properties.VariableNames);
                sel_vars = startsWith(vars, feature_names);
                input_table = input_table(:,sel_vars);
            end
            
            [norm_data, ~, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var);
            
            switch method
                case "UMAP"
%                     n_neighbors = 15;
                    color_file = fullfile(rg.Recordings(1).getParentPath(),'umap','colorsByName.properties');
                    [reduction, umap, clusterIdentifiers, extras] = run_umap(norm_data,'n_components',n_dims,'n_neighbors',10,'min_dist',0.1,'cluster_detail','adaptive','spread',1,'sgd_tasks',20,...
                        'verbose','none','color_file',color_file);
                    rg.DimensionalityReduction.(level).(method).Graph = umap.search_graph;
%                     [M,Q] = community_louvain(umap.search_graph, 0.2);
%                     t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all','Reproducible',true);
%                     clf = fitcensemble(norm_data, M, 'Method','Bag','NumLearningCycles',500,'Learners',t,'Options',statset("UseParallel",true));
%                     explainer = shapley(clf,'QueryPoint',norm_data(1,:),'UseParallel',true);
%                     shap = fit(explainer,norm_data(1,:),'UseParallel',true,'MaxNumSubsets',1000);
%                     [reduction, umap, clusterIdentifiers, extras] = run_umap(norm_data,'n_components',n_dims,'sgd_tasks',20,...
%                         'verbose','none','color_file',color_file);
                case "PCA"
                    [coeff,reduction,latent,tsquared,explained,mu] = pca(norm_data);
                    % abs_coefs = abs(coeff(:,1));
                    % [sorted, sort_idx] = sort(abs_coefs,'descend');
                    % N_feat = 5;
                    % figure; biplot(coeff(sort_idx(1:N_feat),1:2),'Scores',reduction(:,1:2),'VarLabels',feature_names(sort_idx(1:N_feat)))
                case "tSNE"
                    reduction = tsne(norm_data,'Algorithm','exact');%,'Exaggeration',10,'Perplexity',30);
            end
            reduction = reduction(:,1:n_dims);
            
            rg.DimensionalityReduction.(level).(method).Reduction = reduction;
            rg.DimensionalityReduction.(level).(method).GroupingVariable = grouping_var;
            rg.DimensionalityReduction.(level).(method).UnitFeatures = unit_features;
            rg.DimensionalityReduction.(level).(method).NetworkFeatures = network_features;
            rg.DimensionalityReduction.(level).ObjectGroup = object_group;
%             rg.plot_dimensionality_reduction(reduction);
        end
        
        function result = predictAge(rg, level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, ...
                useClustered, normalization_var, N_hyper, K_fold)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                alg string = "rf" %rf, svm, cnb, knn
                stratification_var = "Mutation" %Specify the variable by which to evenly split training and test dataset
                stratification_values = [] %Corresponding to stratification_var, if a value is specified then this will be used as training data (e.g. train on wildtype, test on mutation)
                pooling_vals cell = {}
                network_features string = "all"
                unit_features string = "all"
                useClustered logical = false
                normalization_var string = "PlatingDate"
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
            end
            
            if any(arrayfun(@(x) ~isfield(rg.Recordings(x).Metadata,"DIV"),1:length(rg.Recordings)))
                error('Age data missing')
            end
            
            object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);
            
            if isempty(stratification_values)
                stratification_var = [stratification_var "DIV"];
                [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var);
                if K_fold == -1 %set kfold to LOOCV
                    K_fold = length(rec_group_idx);
                end
                cv = cvpartition(rec_group_idx,'KFold',K_fold);
                t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');
                
                for k = 1:K_fold
                    if level == "Unit"
                        train_table = object_group(cv.training(k)).getUnitFeatures(unit_features);
                        test_table = object_group(cv.test(k)).getUnitFeatures(unit_features);
                        train_idx = logical([ones(1, size(train_table,1)) zeros(1, size(test_table,1))]);
                        test_idx = ~train_idx;
                        input_group = [object_group(cv.training(k)).Units, object_group(cv.test(k)).Units];
                        input_table = [train_table;test_table];
                        unit_recordings = [object_group.MEArecording];
                        metadata = [unit_recordings.Metadata];
                        true_age = [metadata.DIV];
                        Y_train = true_age(train_idx);
                        Y_test = true_age(test_idx);
                    else
                        input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                        train_idx = cv.training(k);
                        test_idx = cv.test(k);
                        input_group = object_group;
                        metadata = [object_group.Metadata];
                        true_age = [metadata.DIV];
                        Y_train = true_age(train_idx);
                        Y_test = true_age(test_idx);
                    end
                    
                    [X_train, X_test, feature_names] = prepareInputMatrix(rg, input_table, input_group, normalization_var, train_idx, test_idx);
                    
                    clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t);
                    Y_pred = predict(clf, X_test);
                    
                    options = statset('UseParallel',true);
                    predImp = clf.oobPermutedPredictorImportance('Options',options);
                    %                     predImp = [];
                    
                    result(k).Mdl = clf;
                    result(k).Y_pred = Y_pred;
                    result(k).Y_test = Y_test';
                    result(k).mse_train = resubLoss(clf);
                    result(k).objects = object_group(test_idx);
                    result(k).predImp = predImp;
                    result(k).GroupLabels = group_labels_comb;
                    result(k).feature_names = feature_names;
                end
            else
                train_group = find(group_labels_comb == stratification_values);
                
                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [object_group.Units];
                    [unit_group_idx, ~] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);
                    train_idx = unit_group_idx == train_group;
                    test_idx = unit_group_idx ~= train_group;
                    unit_recordings = [object_group.MEArecording];
                    metadata = [unit_recordings.Metadata];
                    true_age = [metadata.DIV];
                    
                elseif level == "Recording"
                    input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                    train_idx = rec_group_idx == train_group;
                    test_idx = rec_group_idx ~= train_group;
                    metadata = [object_group.Metadata];
                    true_age = [metadata.DIV];
                    
                end
                [X_train, X_test] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);
                Y_train = true_age(train_idx);
                Y_test = true_age(test_idx);
                t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');
                clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t,'Options',statset('UseParallel',true));
                Y_pred = predict(clf,X_test);
                
                options = statset('UseParallel',true);
                predImp = clf.oobPermutedPredictorImportance('Options',options);
%                 predImp = [];
                
                result.Mdl = clf;
                result.Y_pred = Y_pred;
                result.Y_test = Y_test';
                result.mse_train = resubLoss(clf);
                result.objects = object_group(test_idx);
                result.predImp = predImp;
                result.GroupLabels = group_labels_comb;
            end
            
            rg.Regression.DIV = result;
        end
        
        function result = classifyByFeatureGroups(rg, level, alg, classification_var, pooling_vals, network_features, unit_features, feature_names, useClustered, ... 
                                    grouping_var, grouping_values, normalization, normalization_var, N_hyper, K_fold, tolerance)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                alg string = "rf" %rf, svm, 
                classification_var string = "Mutation" %Metadata field that determines y_test/y_train
                pooling_vals cell = {} %Values corresponding to classification_var that is being classified for(e.g. "LNA")
                network_features string = "all"
                unit_features string = "all"
                feature_names string = []
                useClustered logical = false
                grouping_var string = [] %Metadata field that groups recordings to cultures
                grouping_values = nan %Selected values corresponding to grouping_var
                normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
                normalization_var string = []%"PlatingDate" % normalization by each value of normalization_var
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = -1 % number of K-fold CV
                tolerance double = 1
            end
            
            if isempty(grouping_var) %Get group indices for stratification
                object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            end
            [group_idx, group_labels_comb] = rg.combineMetadataIndices(object_group, classification_var, pooling_vals);
            
            
            if K_fold == -1 %set kfold to LOOCV
                K_fold = length(group_idx);
            end
            cv = cvpartition(group_idx,'KFold',K_fold);
            
            if isempty(grouping_var)
                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [rg.Units];
                    [Y, group_labels_comb] = combineMetadataIndices(rg, object_group, classification_var, pooling_vals);
                    if ~isempty(classification_val) %New division / pool conditions
                        [Y,group_labels_comb] = poolMetadataValues(clf_group_idx, group_labels_comb, classification_val);
                    end
                elseif level == "Recording"
                    input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                    Y = group_idx;
                else 
                    error('Unknown level')
                end
            else
                if level == "Unit"
                    Y = arrayfun(@(x) ones(1,length([object_group{x}(1).Units]))*group_idx(x),1:length(group_idx),'un',0);
                    object_group = cellfun(@(x) [x(1).Units],object_group,'un',0);
                    object_group = [object_group{:}];
                else
                    Y = group_idx;
                end
            end
            
            if ~isempty(feature_names)
                vars = string(input_table.Properties.VariableNames);
                sel_vars = startsWith(vars, feature_names);
                input_table = input_table(:,sel_vars);
            end
            
            for k = 1:K_fold
                [Y_train, Y_test, train_idx, test_idx] = rg.cv_split(Y, cv, k);
                
                [X_train, X_test, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);
                
                [clf,train_acc] = rg.create_classifier(X_train, Y_train, alg, N_hyper);
                [Y_pred,scores] = predict(clf, X_test);
                
                options = statset('UseParallel',true);
                predImp = oobPermutedPredictorImportance(clf,'Options',options);
% predImp = [];
                
                result(k).Mdl = clf;
                result(k).Y_pred = Y_pred;
                result(k).Y_test = Y_test;
                result(k).scores = scores;
                result(k).objects = object_group(test_idx);
                result(k).train_acc = train_acc;
                result(k).predImp = predImp;
                result(k).GroupLabels = group_labels_comb;
                result(k).Features = feature_names;
            end
            rg.Classification.(classification_var) = result;
        end
                
        function [train_acc, test_acc, avg_score, avg_pred_imp] = assessClassifier(rg, clf_result)
            arguments
               rg RecordingGroup
               clf_result struct
            end
            y_pred = vertcat(clf_result.Y_pred);
            y_test = horzcat(clf_result.Y_test)';
            scores = vertcat(clf_result.scores);
            pred_imp = vertcat(clf_result.predImp);
                        
            idx = sub2ind(size(scores), 1:size(scores,1),y_test');
            
            test_acc = sum(y_pred == y_test)/length(y_pred);
            train_acc = mean([clf_result.train_acc]);
            avg_score = mean(scores(idx));
            avg_pred_imp = mean(pred_imp);
            
            fprintf('Training accuracy was %.2f\n', train_acc)
            fprintf('Test accuracy was %.2f\n', test_acc)
            fprintf('Average score was %.2f\n', avg_score)
        end
        
        function result = applyClassifier(rg, level, alg, classification_var, clf_pooling_vals, test_var, test_pooling_vals, network_features, unit_features, feature_names, useClustered,...
                                    grouping_var, grouping_values, normalization, normalization_var, N_hyper, tolerance)
           arguments
               rg RecordingGroup
               level string = "Recording" %Unit or Recording
               alg string = "rf" %rf, svm,
               classification_var string = "Mutation" %Metadata field that determines y_test/y_train
               clf_pooling_vals cell = {}
               test_var string = "Treatment"
               test_pooling_vals cell = {} %First cell corresponds to training data, second to test data
               network_features string = "all"
               unit_features string = "all"
               feature_names string = []
               useClustered logical = false
               grouping_var string = [] %Metadata field that groups recordings to cultures
               grouping_values = nan %Selected values corresponding to grouping_var
               normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
               normalization_var string = "PlatingDate" % normalization by each value of normalization_var
               N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
               tolerance double = 1
           end
           if isempty(grouping_var)
               if level == "Unit"
                   object_group = [rg.Units];
                   input_table = object_group.getUnitFeatures(unit_features);
               elseif level == "Recording"
                   object_group = [rg.Recordings];
                   input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
               else
                   error('Unknown level')
               end
           else
               [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
               if level == "Unit"
                   object_group = cellfun(@(x) [x(1).Units],object_group,'un',0);
                   object_group = [object_group{:}];
               elseif level == "Recording"
                   [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
               else
                   error('Unknown level')
               end
           end
           
           if ~isempty(feature_names)
               vars = string(input_table.Properties.VariableNames);
               sel_vars = startsWith(vars, feature_names);
               input_table = input_table(:,sel_vars);
           end
            
           % Split by classification variable
            [Y, group_labels_comb] = combineMetadataIndices(rg, object_group, classification_var, clf_pooling_vals);
            
            % Identify training and test set
            [test_group_idx, test_group_labels] = combineMetadataIndices(rg, object_group, test_var, test_pooling_vals);
%             [test_group_idx,test_group_labels] = rg.poolMetadataValues(test_group_idx, test_group_labels, );
            
            train_idx = test_group_idx == 1;
            test_idx = ~train_idx;
            Y_train = Y(train_idx);
            Y_test = Y(test_idx);
            
            [X_train, X_test] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);
            [clf,train_acc] = rg.create_classifier(X_train, Y_train, alg, N_hyper);
            [Y_pred,scores] = predict(clf, X_test);
            result.Mdl = clf;
            result.Y_pred = Y_pred;
            result.Y_test = Y_test;
            result.scores = scores;
            result.objects = object_group(test_idx);
            result.train_acc = train_acc;
            result.GroupLabels = group_labels_comb;
            
            rg.Classification.(classification_var).(test_var) = result;
        end
        
        function accuracy_mat = assessAppliedClassifier(rg, result, assessment_var, pooling_vals)
            arguments
                rg RecordingGroup
                result struct %Result of applyClassifier function
                assessment_var string %Metadata variable (e.g. "Treatment")
                pooling_vals cell = {} 
            end
            
            [Y_assess, group_labels_comb] = rg.combineMetadataIndices(result.objects, assessment_var, pooling_vals);
            
            accuracy_mat = nan(length(result.GroupLabels),length(group_labels_comb));
            for t = 1:length(result.GroupLabels)
                for a = 1:length(group_labels_comb)
                    comp_idx = (Y_assess == a) & (result.Y_test == t);
                    accuracy_mat(t,a) = sum(result.Y_test(comp_idx) == result.Y_pred(comp_idx)) ./ sum(comp_idx);
                end
            end
            figure('Color','w');
            heatmap(accuracy_mat,'ColorLimits',[0 1],'XDisplayLabels',group_labels_comb,'YDisplayLabels',result.GroupLabels)
        end
        
        function [train_acc, test_acc, avg_score] = classifyByFeatureGroupsAndGroupingVar(rg, level, alg, classification_var, pooling_vals, useClustered,...
                grouping_var, grouping_value_list, normalization, normalization_var, N_hyper, K_fold, tolerance)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                alg string = "rf" %rf, svm,
                classification_var string = "Mutation" %Metadata field that determines y_test/y_train
                pooling_vals cell = {} %Values corresponding to classification_var that is being classified for(e.g. "LNA")
                useClustered logical = false
                grouping_var string = "DIV" %Metadata field that groups recordings to cultures
                grouping_value_list = nan %Selected values corresponding to grouping_var
                normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
                normalization_var string = "PlatingDate" % normalization by each value of normalization_var
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
                tolerance double = 1
            end
            % Find all existing feature names
             all_unit_features = {"ActivityFeatures", "WaveformFeatures", ["RegularityFeatures","Catch22"], "all"};%[string(fieldnames([rg.Recordings.NetworkFeatures])); "all"];
             all_network_features = {["Regularity", "Catch22"], "Burst", "GraphFeatures", "all"};%[string(fieldnames([rg.Recordings.UnitFeatures])); "all"];
            
            train_acc = nan(length(all_network_features) + length(all_unit_features) + 1, length(grouping_value_list));
            test_acc = nan(length(all_network_features) + length(all_unit_features) + 1, length(grouping_value_list));
            avg_score = nan(length(all_network_features) + length(all_unit_features) + 1, length(grouping_value_list));
            
            % First iterate over network features
            unit_features = [];
            for iNW = 1:length(all_network_features)
                for iGV = 1:length(grouping_value_list)
                    result = classifyByFeatureGroups(rg, level, alg, classification_var, pooling_vals, all_network_features{iNW}, unit_features, [], useClustered,...
                        grouping_var, grouping_value_list(iGV), normalization, normalization_var, N_hyper, K_fold, tolerance);
                    [train_acc(iNW, iGV), test_acc(iNW, iGV), avg_score(iNW, iGV)] = assessClassifier(rg, result);
                end
            end
            
            % Then over all unit features
            network_features = [];
            for iU = 1:length(all_unit_features)
                for iGV = 1:length(grouping_value_list)
                    result = classifyByFeatureGroups(rg, level, alg, classification_var, pooling_vals, network_features, all_unit_features{iU}, [], useClustered,...
                        grouping_var, grouping_value_list(iGV), normalization, normalization_var, N_hyper, K_fold, tolerance);
                    [train_acc(iU + iNW, iGV), test_acc(iU + iNW, iGV), avg_score(iU + iNW, iGV)] = assessClassifier(rg, result);
                end
            end
            
            % Finally over all feature groups combined
            network_features = "all";
            unit_features = "all";
            for iGV = 1:length(grouping_value_list)
                result = classifyByFeatureGroups(rg, level, alg, classification_var, pooling_vals, network_features, unit_features, [], useClustered,...
                    grouping_var, grouping_value_list(iGV), normalization, normalization_var, N_hyper, K_fold, tolerance);
                [train_acc(end, iGV), test_acc(end, iGV), avg_score(end, iGV)] = assessClassifier(rg, result);
            end
            rg.Classification.FeatureGroupsByGroupingVar.TrainAcc = train_acc;
            rg.Classification.FeatureGroupsByGroupingVar.TestAcc = test_acc;
            rg.Classification.FeatureGroupsByGroupingVar.AvgScore = avg_score;
            rg.Classification.FeatureGroupsByGroupingVar.FeatureGroups = [all_network_features all_unit_features {"combined"}];
            rg.Classification.FeatureGroupsByGroupingVar.GroupingVar = grouping_var;
            rg.Classification.FeatureGroupsByGroupingVar.GroupingValueList = grouping_value_list;
        end
        
        function result = regressionByFeatureGroups(rg, level, regression_var, stratification_var, stratification_values, pooling_vals, grouping_var, grouping_values, ...
                network_features, unit_features, useClustered, normalization_var, normalization, tolerance, N_hyper, K_fold)
            arguments
                rg RecordingGroup
                level string = "Recording" %Unit or Recording
                regression_var string = "Concentration" %Metadata variable that will be regressed for
                stratification_var = "Mutation" %Specify the variable by which to split training and test dataset (e.g. train on wildtype, test on mutation)
                stratification_values = [] %Corresponding to stratification_var, if a value is specified then this will be used as training data
                pooling_vals cell = {}
                grouping_var = "DIV"
                grouping_values = nan
                network_features string = "all"
                unit_features string = "all"
                useClustered logical = false
                normalization_var string = "PlatingDate"
                normalization string = []
                tolerance double = 1
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
            end
            
            if any(arrayfun(@(x) ~isfield(rg.Recordings(x).Metadata,regression_var),1:length(rg.Recordings)))
                error('Regression data missing')
            end
            
            [~, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);
            
            if isempty(stratification_values)
                stratification_var = [stratification_var regression_var];
                [rec_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);
                if K_fold == -1 %set kfold to LOOCV
                    K_fold = length(rec_group_idx);
                end
                cv = cvpartition(rec_group_idx,'KFold',K_fold);
                t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');
                
                for k = 1:K_fold
                    if level == "Unit"
                        error('Not yet implemented')
%                         train_table = object_group(cv.training(k)).getUnitFeatures(unit_features);
%                         test_table = object_group(cv.test(k)).getUnitFeatures(unit_features);
%                         train_idx = logical([ones(1, size(train_table,1)) zeros(1, size(test_table,1))]);
%                         test_idx = ~train_idx;
%                         input_group = [object_group(cv.training(k)).Units, object_group(cv.test(k)).Units];
%                         input_table = [train_table;test_table];
%                         unit_recordings = [object_group.MEArecording];
%                         metadata = [unit_recordings.Metadata];
%                         true_value = [metadata.(regression_var)];
%                         Y_train = true_value(train_idx);
%                         Y_test = true_value(test_idx);
                    else
                        %                         input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                        [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
                        train_idx = cv.training(k);
                        test_idx = cv.test(k);
                        [reg_group_idx, group_labels_comb] = combineMetadataIndices(rg, object_group, regression_var, pooling_vals);
                        true_value = group_labels_comb(reg_group_idx);
                        if isstring(true_value)
                           true_value = str2double(true_value); 
                        end
                        input_group = object_group;
                        Y_train = true_value(train_idx);
                        Y_test = true_value(test_idx);
                    end
                    
                    [X_train, X_test] = prepareInputMatrix(rg, input_table, input_group, normalization_var, train_idx, test_idx);
                    
                    clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t,'Options',statset('UseParallel',true));
                    Y_pred = predict(clf, X_test);
                    %                     predImp = clf.oobPermutedPredictorImportance;
                    predImp = [];
                    
                    result(k).Mdl = clf;
                    result(k).Y_pred = Y_pred;
                    result(k).Y_test = Y_test;
                    result(k).mse_train = resubLoss(clf);
                    result(k).objects = object_group(test_idx);
                    result(k).predImp = predImp;
                    result(k).GroupLabels = group_labels_comb;
                end
            else
                error('Not yet implemented')
                train_group = find(group_labels_comb == stratification_values);
                
                if level == "Unit"
                    input_table = object_group.getUnitFeatures(unit_features);
                    object_group = [object_group.Units];
                    [unit_group_idx, ~] = combineMetadataIndices(rg, object_group, stratification_var, pooling_vals);
                    train_idx = unit_group_idx == train_group;
                    test_idx = unit_group_idx ~= train_group;
                    unit_recordings = [object_group.MEArecording];
                    metadata = [unit_recordings.Metadata];
                    true_value = [metadata.(regression_var)];
                    
                elseif level == "Recording"
                    input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                    train_idx = rec_group_idx == train_group;
                    test_idx = rec_group_idx ~= train_group;
                    metadata = [object_group.Metadata];
                    true_value = [metadata.(regression_var)];
                    
                end
                [X_train, X_test] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);
                Y_train = true_value(train_idx);
                Y_test = true_value(test_idx);
                t = templateTree('Surrogate','on','MinLeafSize',1,'NumVariablesToSample','all');
                clf = fitrensemble(X_train, Y_train,'Method','Bag','NumLearningCycles',500,'Learners',t,'Options',statset('UseParallel',true));
                Y_pred = predict(clf,X_test);
                %                     predImp = clf.oobPermutedPredictorImportance;
                predImp = [];
                
                result(k).Mdl = clf;
                result(k).Y_pred = Y_pred;
                result(k).Y_test = Y_test;
                result(k).mse_train = resubLoss(clf);
                result(k).objects = object_group(test_idx);
                result(k).predImp = predImp;
                result(k).GroupLabels = group_labels_comb;
                
            end
            rg.Regression.(regression_var) = result;
        end
        
        function [iMetadata, metadata_groups] = returnMetadataArray(rg, metadata_object, metadata_name)
           arguments
               rg RecordingGroup
               metadata_object %Can be a string referring to the respective group (Unit,Recording,Culture) within the RecordingGroup, or a Unit,Recording or Culture array
               metadata_name string
           end
           
           if isstring(metadata_object)
              switch metadata_object
                  case "Unit"
                      metadata_object = [rg.Units.MEArecording];
                      
                  case "Recording"
                      metadata_object = [rg.Recordings];
                      
                  case "Culture"
                      metadata_object = cellfun(@(x) x(1),rg.Cultures);
                      
              end
              
           elseif iscell(metadata_object)
               metadata_object = cellfun(@(x) x(1),metadata_object);
               
           elseif class(metadata_object) == "Unit"
               metadata_object = [metadata_object.MEArecording];
           end
           
           metadata_struct = [metadata_object.Metadata];
           assert(isfield(metadata_struct,metadata_name),"Not a valid metadata field")
           metadata = [metadata_struct.(metadata_name)];
           [iMetadata, metadata_groups] = findgroups(metadata);
        end
        
        function [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, metadata_names, pooling_vals)
            arguments
                rg RecordingGroup
                metadata_object %Can be a string referring to the respective group (Unit,Recording,Culture) within the RecordingGroup, or a Unit,Recording or Culture array
                metadata_names string
                pooling_vals cell = {}
            end
            
            idx_array = cell(1,length(metadata_names));
            group_array = cell(1,length(metadata_names));
            for m = 1:length(metadata_names)
                [idx, group] = returnMetadataArray(rg, metadata_object, metadata_names(m)); %{end-m+1}
                if ~isempty(pooling_vals)
                    new_idx_array = [];
                    new_groups = string();
                    for p = 1:length(pooling_vals{m})
                        pool_idx = find(matches(string(group), string(pooling_vals{m}{p})));
                        new_idx_array(p,:) = matches(string(idx), string(pool_idx)) * p;
                        new_groups(p) = join(string(group(pool_idx)),' ');
                    end
                    idx_array{end-m+1} = sum(new_idx_array);
                    group_array{end-m+1} = new_groups;
                else
                    idx_array{end-m+1} = idx;
                    group_array{end-m+1} = group;
                end
            end
            group_vals = cellfun(@unique, idx_array, "un",0);
            group_combs = combvec(group_vals{:});
            group_idx = vertcat(idx_array{:});
            [~,cluster_idx] = ismember(group_idx', group_combs','rows');
            group_labels = arrayfun(@(x) group_array{x}(group_combs(x,:))',1:length(group_array),"un",0);
            if length(metadata_names) == 1
                group_labels_comb = [group_labels{:}];
            else
                group_labels_comb = join(string([group_labels{end:-1:1}]));
                cluster_ids = unique(cluster_idx);
                max_clust = length(group_labels_comb);
                group_labels_comb = group_labels_comb(cluster_ids);
                new_ids = cumsum(sum(cluster_ids == 1:max_clust)); %Remove ids without corresponding group
                cluster_idx = new_ids(cluster_idx);
            end
            if ~isstring(group_labels_comb)
                group_labels_comb = string(group_labels_comb);
            end
        end
        
        function culture = findCulture(rg, rec)
            arguments
                rg RecordingGroup
                rec MEArecording
            end
            culture_idx = cellfun(@(x) any(ismember(x, rec)),rg.Cultures);
            culture = rg.Cultures{culture_idx};
        end

        function [sel_recs, group_labels_comb] = findBestRecordings(rg, metadata_names, pooling_vals)
            arguments
                rg RecordingGroup
                metadata_names string
                pooling_vals cell = {}
            end
            if ~isempty(metadata_names)
                [cluster_idx, group_labels_comb] = rg.combineMetadataIndices(rg.Recordings, metadata_names, pooling_vals);
            else
                cluster_idx = ones(1,length(rg.Recordings));
                group_labels_comb = [];
            end
            N_units = cellfun(@length, {rg.Recordings.Units});

            sel_recs = MEArecording();
            for c = 1:length(unique(cluster_idx))
                sel_idx = find(cluster_idx == c);
                max_idx = sel_idx(N_units(sel_idx) == max(N_units(sel_idx)));
                sel_recs(c) = rg.Recordings(max_idx(1));
            end
        end

        function feature_names = returnFeatureNames(rg, feature_groups, feature_subgroups)
            arguments
               rg RecordingGroup
               feature_groups string = ["UnitFeatures","NetworkFeatures"]
               feature_subgroups cell  = {}
            end
            
            
            for fg = 1:length(feature_groups)
                if isempty(feature_subgroups)
                    feature_sg = string(fieldnames([rg.Recordings.(feature_groups(fg))]));
                else
                    feature_sg = feature_subgroups{fg};
                end
                features = [];
                for fs = 1:length(feature_sg)
                    features = [features string(rg.Recordings(1).(feature_groups(fg)).(feature_sg(fs)).Properties.VariableNames)];
                end
                feature_names{fg} = features;
            end
        end
    
        function [cluster_idx, cmap] = clusterByFeatures(rg, input_type, level, method, N_clust, grouping_var)
           arguments
                rg RecordingGroup
                input_type string = "UMAP" %"RAW","UMAP","PCA"
                level string = "Unit" %Unit, Recording, Culture
                method string = "kmeans" %kmeans, hierarchical, spectral, gmm, louvain
                N_clust {isnumeric} = nan %Number of imposed clusters // nan selects N_clust by finding the optimal silhouette score
                grouping_var string = "PlatingDate" % only for input_type "RAW"
           end
           
           if strcmp(input_type, "RAW")
               object_group = [rg.Units];
               unit_feature_array = [object_group.AverageFeatures];
               feature_table = rg.Recordings.concatenateFeatures(unit_feature_array, "all");
               input_mat = feature_table.Variables;
               input_mat(isnan(input_mat)) = 0;%Handle rare case where NaN appears
               norm_data = normalizeByGroup(rg, input_mat, object_group, grouping_var); %Normalize
               feat_mat = norm_data./max(abs(norm_data)); %Scale data
           elseif contains(input_type, ["PCA","UMAP"])
               feat_mat = rg.DimensionalityReduction.(level).(input_type).Reduction;
           else
              error("Unknown input type") 
           end
           
           switch method
               case "kmeans"
                   cluster_fun = @(x,k) kmeans(x,k); 
                   
               case "hierarchical"
                   cluster_fun = @(x,k) clusterdata(x,k); 
                   
               case "dbscan"
                   error("Not yet implemented")
                   
               case "spectral"
                   cluster_fun = @(x,k) spectralcluster(x,k);
                   
               case "gmm"
                   assert(size(feat_mat,1)>size(feat_mat,2),"GMM needs more samples than variables")
                   cluster_fun = @(x,k) cluster(fitgmdist(x,k),x);
                   
               case "louvain"
                   assert(input_type == "UMAP","UMAP graph needed for louvain clustering")
                   [M,Q] = community_louvain(rg.DimensionalityReduction.(level).UMAP.Graph, N_clust); %N_Clust refers here to the resolution parameter
                   
           end
           if isnan(N_clust)
               evaluation = evalclusters(feat_mat,cluster_fun,'silhouette',"KList",1:10);
               N_clust = evaluation.OptimalK;
               cluster_idx = cluster_fun(feat_mat,N_clust);
           elseif method == "louvain"
               cluster_idx = M;
%                fprintf('Modularity: %.2f\n',Q)
           else
               cluster_idx = cluster_fun(feat_mat,N_clust);
           end 
           rg.Clustering.(level).(method).Index = cluster_idx;
           rg.Clustering.(level).(method).k = N_clust;
           rg.Clustering.(level).(method).Input = input_type;
           
%            figure('Color','w');
           [reduction, cmap] = RecordingGroup.plot_cluster_outlines(feat_mat, cluster_idx, gca);
           title(sprintf('%s clustering on %ss',method, level))
        end
        
        function cluster_purity = calculateClusterPurity(rg, input_type, method, cluster_vars, pooling_vals)
           arguments
               rg RecordingGroup
               input_type string = "UMAP"
               method string = "kmeans" %Clustering method
               cluster_vars string = "Mutation" %Metadata variable for the ground truth
               pooling_vals cell = {}
           end
           
           [ground_truth_labels, ~] = combineMetadataIndices(rg, "Culture", cluster_vars, pooling_vals);
           N_clust = max(ground_truth_labels);
           if method == "louvain"
               predicted_labels = 0;
               res_param = 0.9;
               while max(predicted_labels) < N_clust
                   predicted_labels = rg.clusterByFeatures(input_type, "Recording", method, res_param);
                   res_param = res_param + 0.01;
                   if max(predicted_labels) > N_clust
                       fprintf('Resolution parameter: %.2f ',res_param)
                       error('Could not find suitable resolution parameter')
                   end
               end
           else
               predicted_labels = rg.clusterByFeatures(input_type, "Recording", method, N_clust);
           end
           cm = confusionmat(ground_truth_labels,predicted_labels);

           cluster_purity = sum(max(cm)) / length(predicted_labels);
        end
        
        function result = classifyClusteredUnits(rg, true_cluster_idx, feature_groups, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance)
            arguments
                rg RecordingGroup
                true_cluster_idx double
                feature_groups string
                grouping_var string = [] %Metadata field that groups recordings to cultures
                grouping_values = nan %Selected values corresponding to grouping_var
                feature_names string = []
                normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
                normalization_var string = "PlatingDate" % normalization by each value of normalization_var
                N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
                K_fold (1,1) double = 5 % number of K-fold CV
                tolerance double = 1
            end
            network_features = [];
            useClustered = false;
            level = "Unit";
            alg = "rf";
            
            if isempty(grouping_var)
                object_group = [rg.Recordings];
                input_table = object_group.getUnitFeatures(feature_groups);
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, feature_groups, normalization, useClustered);
           
            end
            
            cv = cvpartition(true_cluster_idx,'KFold',K_fold);
            
            if ~isempty(feature_names)
                vars = string(input_table.Properties.VariableNames);
                sel_vars = startsWith(vars, feature_names);
                input_table = input_table(:,sel_vars);
            end
            
            for k = 1:K_fold
                [Y_train, Y_test, train_idx, test_idx] = rg.cv_split(true_cluster_idx, cv, k);
                
                % Balance class distribution by oversampling
                if oversample
                    [input_table, Y_train, train_idx, test_idx] = oversample_smaller_class(input_table, Y_train, train_idx);
                end
                %
                [X_train, X_test, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);
                
                [clf,train_acc] = rg.create_classifier(X_train, Y_train, alg, N_hyper);
                [Y_pred,scores] = predict(clf, X_test);
                
                options = statset('UseParallel',true);
                predImp = oobPermutedPredictorImportance(clf,'Options',options);
% predImp = [];
                
                result(k).Mdl = clf;
                result(k).Y_pred = Y_pred;
                result(k).Y_test = Y_test;
                result(k).scores = scores;
                result(k).train_acc = train_acc;
                result(k).predImp = predImp;
                result(k).Features = feature_names;
            end
            rg.Classification.UnitClusters = result;
        end
        
        function result = applyClusteredUnitClassifier(rg, true_cluster_idx, train_idx, test_idx, feature_groups, clf,grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance)
            arguments
               rg RecordingGroup
               true_cluster_idx double %Only needed if clf is not provided
               train_idx logical
               test_idx logical
               feature_groups string
               clf = []
               grouping_var string = [] %Metadata field that groups recordings to cultures
               grouping_values = nan %Selected values corresponding to grouping_var
               feature_names string = []
               normalization = [] %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
               normalization_var string = "PlatingDate" % normalization by each value of normalization_var
               N_hyper (1,1) double = 0 %If >0 do hyperparameter optimization
               tolerance double = 0
            end
            network_features = [];
            useClustered = false;
            level = "Unit";
            alg = "rf";
            
            if isempty(grouping_var)
                object_group = [rg.Recordings];
                input_table = object_group.getUnitFeatures(feature_groups);
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, feature_groups, normalization, useClustered);
           
            end
            
            if ~isempty(feature_names)
                vars = string(input_table.Properties.VariableNames);
                sel_vars = matches(vars, feature_names);
                input_table = input_table(:,sel_vars);
            end
            
            [X_train, X_test, feature_names] = prepareInputMatrix(rg, input_table, object_group, normalization_var, train_idx, test_idx);
            
            if isempty(clf)
                [clf,train_acc] = rg.create_classifier(X_train, true_cluster_idx, alg, N_hyper);
                result.train_acc = train_acc;
            end
            
            [Y_pred,scores] = predict(clf, X_test);
            
            result.Y_pred = Y_pred;
            result.scores = scores;
            result.feature_names = feature_names;
            result.Mdl = clf;
        end
        
        function assignUnitClusterIdx(rg,method,calc_feat, concatenated)
           arguments
               rg RecordingGroup
               method string = "louvain"
               calc_feat logical = true %(Re)calculate unit feature averages per cluster
               concatenated logical = false %Check if recordings were concatenated (assigns unit IDs to all concatenated timepoints)
           end
           cluster_idx = rg.Clustering.Unit.(method).Index;
           cluster_idx = num2cell(cluster_idx); %Prepare to use deal to assign cluster ids
           [rg.DimensionalityReduction.Unit.ObjectGroup.ClusterID] = deal(cluster_idx{:});
           
           if calc_feat
               N_clust = num2cell(ones(size(rg.Recordings))*max([cluster_idx{:}]));
               [rg.Recordings.NumUnitClusters] = deal(N_clust{:});
               if concatenated
                  for c = 1:length(rg.Cultures)
                     ids = num2cell([rg.Cultures{c}(1).Units.ClusterID]); %Prepare to use deal to assign cluster ids
                     for r = 2:length(rg.Cultures{c})
                         [rg.Cultures{c}(r).Units.ClusterID] = deal(ids{:});
                     end
                  end
               end
               for r = 1:length(rg.Recordings)
                  rg.Recordings(r).calculateClusterSingleCellFeatures();
               end
           end
        end
        
        function removeUnitsByCluster(rg, method, ID, concatenated)
            arguments
               rg RecordingGroup
               method string
               ID double
               concatenated = false % If true, units are removed for each recording of a culture and not only for individual recordings
            end
            rg.assignUnitClusterIdx(method,false);
            
            if concatenated
                for c = 1:length([rg.Cultures])
                    rm_idx = [];
                    for r = 1:length([rg.Cultures{c}])
                        rm_idx = [rm_idx find([rg.Cultures{c}(r).Units.ClusterID] == ID)];
                    end
                    rm_idx = unique(rm_idx);
                    for r = 1:length([rg.Cultures{c}])
                        rg.Cultures{c}(r).Units(rm_idx) = [];
                    end
                end
            else
                
                for r = 1:length([rg.Recordings])
                    units = [rg.Recordings(r).Units];
                    rg.Recordings(r).Units([units.ClusterID] == ID) = [];
                end
            end
            rg.Units = [rg.Recordings.Units];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [cluster_idx, group_labels_comb] = plot_true_clusters(rg, level, method, grouping_var, options)
            arguments
                rg RecordingGroup
                level string = "Unit" %Unit, Recording, Culture
                method string = "UMAP" %dim reduction method: "UMAP","PCA"
                grouping_var string = "Mutation"
                options.pooling_vals cell = {}
                options.plot_centroid = false
                options.nodeSz (1,1) {isnumeric} = 10
                options.mapSz {isnumeric} = 300
                options.sigma (1,1) {isnumeric} = 5
                options.cmap (:,3) double = []
            end
            reduction = rg.DimensionalityReduction.(level).(method).Reduction;
            metadata_object = rg.DimensionalityReduction.(level).ObjectGroup;
            [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, grouping_var, options.pooling_vals);
            if isempty(options.cmap)
                options.cmap = othercolor('Mrainbow',max(cluster_idx));%Mdarkrainbow
            end
            ax = gca;
            RecordingGroup.plot_cluster_outlines(reduction, cluster_idx, ax, options.plot_centroid, options.nodeSz, options.mapSz, options.sigma, options.cmap);
            if length(group_labels_comb) == 1
                group_labels_comb = strsplit(group_labels_comb, ' ');
            end
            legend(group_labels_comb,"Location","best","Box","off")
        end
               
        function plot_single_cluster(rg, reduction, cluster_idx)
            arguments
                rg RecordingGroup
                reduction {isnumeric}
                cluster_idx (1,:) {isnumeric} = ones(1,size(reduction,1))
            end
            
            if length(cluster_idx)>100
                sz = 5;
            else
                sz = 20;
            end
            
            N_clust = length(unique(cluster_idx));
            
            nexttile
            if size(reduction,2) == 2
                scatter(reduction(:,1),reduction(:,2),sz,cluster_idx,'filled')
            else
                scatter3(reduction(:,1),reduction(:,2),reduction(:,3),sz,cluster_idx,'filled')
            end
            
            c_map = othercolor('Set19',N_clust);
            colormap(c_map)
            
        end
        
        function plot_dimensionality_reduction(rg,reduction,cluster_idx,group_labels_comb)
            arguments
                rg RecordingGroup
                reduction {isnumeric}
                cluster_idx (1,:) {isnumeric} = ones(1,size(reduction,1))
                group_labels_comb string = "Cluster " + unique(cluster_idx)
            end
            
            if length(cluster_idx)>500
                sz = 5;
            else
                sz = 20;
            end
            
            N_clust = length(unique(cluster_idx));
            c_map = othercolor('Spectral9',N_clust);
            figure('Color','w');
            if size(reduction,2) == 2
                nexttile
                scatter(reduction(:,1),reduction(:,2),sz,cluster_idx,'filled','MarkerEdgeColor','k')
                colormap(c_map)
                cb = colorbar; cb.Title.String = "Cluster Index";
                cb.Ticks = linspace(cb.Limits(1) + (N_clust-1)/(2*N_clust), cb.Limits(2) - (N_clust-1)/(2*N_clust), N_clust); cb.TickLabels = 1:N_clust;
                if N_clust > 1
                    xl = xlim; yl = ylim;
                    for i = 1:N_clust
                        nexttile
                        scatter(reduction(cluster_idx==i,1),reduction(cluster_idx==i,2),sz,c_map(i,:),'filled','MarkerEdgeColor','k')
                        title(group_labels_comb(i))
                        xlim(xl); ylim(yl)
                    end
                end
            else
                nexttile
                scatter3(reduction(:,1),reduction(:,2),reduction(:,3),sz,cluster_idx,'filled','MarkerEdgeColor','k')
                colormap(c_map)
                cb = colorbar; cb.Title.String = "Cluster Index";
                cb.Ticks = linspace(cb.Limits(1) + (N_clust-1)/(2*N_clust), cb.Limits(2) - (N_clust-1)/(2*N_clust), N_clust); cb.TickLabels = 1:N_clust;
                if N_clust > 1
                    xl = xlim; yl = ylim; zl = zlim;
                    for i = 1:N_clust
                        nexttile
                        scatter3(reduction(cluster_idx==i,1),reduction(cluster_idx==i,2),reduction(cluster_idx==i,3),sz,c_map(i,:),'filled','MarkerEdgeColor','k')
                        title(group_labels_comb(i))
                        xlim(xl); ylim(yl);zlim(zl)
                    end
                end
            end
        end
        
        function [colors, avg_wf] = plot_cluster_waveforms(rg,method,cluster_idx, object_group)
            arguments
                rg RecordingGroup
                method string = "kmeans"
                cluster_idx {isnumeric} = []
                object_group = []
            end
            if isempty(cluster_idx)
                cluster_idx = [rg.Clustering.Unit.(method).Index];
            end
            %             aligned_wf = rg.alignWaveforms([rg.Units]);

            if isempty(object_group)
                aligned_wf = rg.alignWaveforms([rg.DimensionalityReduction.Unit.ObjectGroup]);
            else
                aligned_wf = rg.alignWaveforms(object_group);
            end
            N_clust = length(unique(cluster_idx));
            avg_wf = cell(1,N_clust);
            cluster_wfs = cell(1,N_clust);
            cluster_size = histcounts(cluster_idx);
            time_vector = (1:size(aligned_wf,1))/20;
            
            for c = 1:N_clust
                cluster_wfs{c} = aligned_wf(:,cluster_idx == c);
                avg_wf{c} = median(cluster_wfs{c},2);
            end
            
            figure("Color","w");
            ax = nexttile;
            avg_time_vector = linspace(min(time_vector),max(time_vector),max([200,length(time_vector)]));
            smoothed_avg = pchip(time_vector,horzcat(avg_wf{:})',avg_time_vector);
%             p = plot(time_vector,horzcat(avg_wf{:}),"LineWidth",2);
            p = plot(avg_time_vector,smoothed_avg,"LineWidth",2);
            ax.ColorOrder = othercolor('Spectral9',max(cluster_idx));
            
            l = legend("N = " + string(cluster_size)); l.Box = "off";
            box off; xlabel("Time [ms]"); ylim([-1 1])
            for c = 1:N_clust
                nexttile
                plot(time_vector,cluster_wfs{c},"LineWidth",0.1,'Color',[p(c).Color,0.1]);
                hold on
                plot(avg_time_vector,smoothed_avg(c,:),"LineWidth",2,'Color',[0 0 0 0.5])
                box off; axis tight; xlabel("Time [ms]"); ylim([-1 1])
            end
            colors = vertcat(p.Color);
        end
        
        function [feature_table, mean_mat, sd_mat, group_labels_comb] = plot_feature_trajectories(rg, level, grouping_var, grouping_values, network_features, unit_features, normalization,...
                feature_names, comp_var, pooling_vals, useClustered, tolerance, colors)
            arguments
                rg RecordingGroup
                level string
                grouping_var string
                grouping_values (1,:) {isnumeric}
                network_features string
                unit_features string
                normalization string = [] %"baseline" or "scaled"
                feature_names string = []
                comp_var string = "Mutation"
                pooling_vals cell = {}
                useClustered logical = false
                tolerance = 0
                colors double =[]
            end
            [feature_table, culture_array, grouping_values] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            
            if ~isempty(feature_names)
                sel_idx = startsWith(feature_table.Properties.VariableNames, feature_names);
                if isempty(sel_idx)
                    error("Could not find " + feature_names)
                end
                feature_table = feature_table(:,sel_idx);
            end
            
            grouping_values = 1:length(grouping_values);
            if normalization == "baseline"
               % grouping_values =  grouping_values(2:end);
               grouping_values = 1:(length(grouping_values)-1);
            end
            
            
            
            N_features = size(feature_table,2)/numel(grouping_values);
            
            if level == "Recording"
                [group_idx, group_labels_comb] = combineMetadataIndices(rg, culture_array, comp_var, pooling_vals);
                N_groups = length(group_labels_comb);
            elseif level == "Unit"
                %                 group_idx = [rg.Units.ClusterID];
                group_idx = [rg.DimensionalityReduction.Unit.ObjectGroup.ClusterID];
                N_groups = length(unique(group_idx));
                group_labels_comb = "Cluster " + [1:N_groups];
            else
                error("Unknown level")
            end
            
            fontsz = 7;
            
            if isempty(colors)
               colors = othercolor('RdBu4',N_groups);
            end
            
            min_diff = min(diff(grouping_values));
            jitter = linspace(-min_diff/10,min_diff/10,N_groups);%jitter = zeros(1,N_groups); %Change to have jitter
            
            mean_mat = nan(N_features, length(grouping_values), N_groups);
            sd_mat = nan(N_features, length(grouping_values), N_groups);
            features = string();
            
            figure('Color','w','Position',[100 100 1500 1000]);
            tiledlayout('flow','TileSpacing','tight')
            for f = 1:N_features
                if ~isempty(feature_names)
                    sel_idx = find(startsWith(feature_table.Properties.VariableNames, feature_names(f)));
                else
                    sel_idx = f:N_features:size(feature_table,2);
                end
                feature_parts = strsplit(feature_table.Properties.VariableNames{sel_idx(1)},"_");
                features(f) = [feature_parts{1:end-1}];
                feature_matrix = feature_table(:,sel_idx).Variables;
                
                nexttile
                
                for g = 1:N_groups
                    x = grouping_values + jitter(g);
                    
                    data_mat = feature_matrix(group_idx==g,:);
                    data_mat(isoutlier(data_mat,'ThresholdFactor',5)) = nan;

                    
                    if normalization == "baseline"
                        x = (1:(length(grouping_values)+1)) + jitter(g);
                        plot_x = 1:length(x);
                        log_data = log(data_mat);
                        log_data(isinf(log_data)) = nan;
                        y = [0 mean(log_data,1,'omitnan')];
                        y_err = [0 std(log_data,[],1,'omitnan')];
                    else
                        y = mean(data_mat,1,'omitnan');
                        y_err = std(data_mat,[],1,'omitnan');
                        plot_x = grouping_values;
                    end
                    xx = linspace(min(x),max(x),((max(x) - min(x))*5)+1);
                    yy = makima(x,y,xx);
                    
                    plot(xx,yy,'Color',colors(g,:),'HandleVisibility','off')
                    hold on
                    errorbar(plot_x+jitter(g),y,y_err,...
                        'LineWidth',1,'Color',colors(g,:),'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',2,...
                        'MarkerFaceColor',colors(g,:));
                    set(gca,'FontSize',fontsz)
                    marg = get(gca,'ylabel');
                    set(marg,'Margin',3)
                    if normalization == "baseline"
                        mean_mat(f,:,g) = y(2:end);
                        sd_mat(f,:,g) = y_err(2:end);
                    else
                        mean_mat(f,:,g) = y;
                        sd_mat(f,:,g) = y_err;
                    end
                end
                l = legend(group_labels_comb,'Box','off');
                
                title(features(f))
                xlabel(grouping_var)
                box off
            end
            
            %             l.Position = [0.02 0.8 0.1 0.1];
        end
        
        function [color_mat, features] = plot_feature_heatmap(rg, level, grouping_var, grouping_values, network_features, unit_features, feature_names, normalization, comp_var, pooling_vals, useClustered, tolerance, color_lim)
           arguments
                rg RecordingGroup
                level string
                grouping_var string
                grouping_values (1,:) {isnumeric}
                network_features string
                unit_features string
                feature_names string
                normalization string = [] %[], "baseline" or "scaled"
                comp_var = "Mutation"
                pooling_vals cell = {}
                useClustered logical = false
                tolerance = 1
                color_lim double = 3 %Maximum value to cap colormap
           end
           
           [feature_table, culture_array, grouping_values] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
           
           if ~isempty(feature_names)
                vars = string(feature_table.Properties.VariableNames);
                sel_vars = startsWith(vars, feature_names);
                feature_table = feature_table(:,sel_vars);
           end
            
           %We return grouping_values in case we extracted the values by
           %setting it to nan 
           [group_idx, group_labels_comb] = combineMetadataIndices(rg, culture_array, comp_var, pooling_vals);
           features = [];
           feat_mat = feature_table.Variables;
           feat_mat(isinf(feat_mat)) = nan;
           norm_mat = normalize(feat_mat,'range',[0 1]);
           means = arrayfun(@(x) mean(norm_mat(group_idx == x,:),'omitnan'),1:length(unique(group_idx)),'un',0);
           color_vector = log(means{1}./means{2});
           color_mat = reshape(color_vector,[],length(grouping_values));
           for f = 1:length(feature_table.Properties.VariableNames)
               feature_parts = strsplit(feature_table.Properties.VariableNames{f},"_");
               feat = string(strjoin(feature_parts(1:end-1),'_'));
               if isempty(features) || ~matches(feat, features)
                   features = [features feat];
               end
           end
           
           if ~isempty(feature_names) %Sort features to match input
               [~,idx] = ismember(feature_names, features);
               color_mat = color_mat(idx,:);
               features = feature_names;
           end
           
           if normalization == "baseline"
              color_mat(:,1) = [];
              grouping_values(1) = [];
           end
           
           figure('Color','w');
           imagesc(color_mat)
           cm = othercolor('RdBu9',100);
           colormap(cm([1:45,55:end],:)) %Make smaller differences more visible
           clim = [log(1/color_lim) log(color_lim)];
           set(gca,'CLim',clim) %Set color limits to 1/3 and 3 
           xticks(1:size(color_mat,2))
           xticklabels(grouping_values)
           xlabel(grouping_var)
           yticks(1:size(color_mat,1))
           yticklabels(features)
           cb = colorbar;
           cb.Location = 'northoutside';
           cb.Title.String = group_labels_comb(1) + "/" + group_labels_comb(2);
           cb.Title.FontWeight = 'bold';
           cb.Ticks = [clim(1) 0 clim(end)];
           cb.TickLabels = [sprintf("<%.1f",1/color_lim) + color_lim, "1", ">" + color_lim];
           set(gca,'FontSize',7)
           set(gca,'TickLength',[0 0])
        end
        
        function [value_array, group_labels_comb] = plot_cluster_densities(rg, level, method, grouping_var, pooling_vals, n_bins, smoothing_factor)
            arguments
                rg RecordingGroup
                level string = "Unit" %Unit, Recording, Culture
                method string = "UMAP" %dim reduction method: "UMAP","PCA"
                grouping_var = "Mutation" %Metadata name or cluster index
                pooling_vals cell = {}
                n_bins (1,1) {isnumeric} = 50
                smoothing_factor double = 0.5
            end
            reduction = rg.DimensionalityReduction.(level).(method).Reduction;
            if isstring(grouping_var)
                metadata_object = rg.DimensionalityReduction.(level).ObjectGroup;
                [cluster_idx, group_labels_comb] = combineMetadataIndices(rg, metadata_object, grouping_var, pooling_vals);
            else
                cluster_idx = grouping_var;
                group_labels_comb = string(unique(cluster_idx));
            end
            N_clust = length(unique(cluster_idx));
            edges = {linspace(min(reduction(:,1)), max(reduction(:,1)),n_bins) linspace(min(reduction(:,2)), max(reduction(:,2)),n_bins)};
            value_array = zeros(n_bins,n_bins,N_clust);
            [X,Y] = meshgrid(edges{1},edges{2});
            figure('Color','w');
            tiledlayout('flow','TileSpacing','compact')
            for i = 1:N_clust
                
                plot_idx = cluster_idx==i;
                data = [reduction(plot_idx,1),reduction(plot_idx,2)];
                ax(i) = nexttile;
                values = hist3(data, edges);
                embedded_values = zeros(size(values)*1.5);
                embed_coor = round(n_bins*0.25)+1:round(n_bins*1.25);
                embedded_values(embed_coor, embed_coor) = values;
                value_array(:,:,i) = values;
                % values = values.'./length(cluster_idx);
% surf(X,Y,values)
                smoothed_data = smoothdata2(embedded_values,"gaussian",SmoothingFactor=smoothing_factor);
                contourf(smoothed_data,'LineWidth',0.25,'FaceAlpha',0.5);
                axis square
                title(group_labels_comb(i))
                xticks([])
                yticks([])
                colormap(flipud(hot))
            end
%             arrayfun(@(x) set(x,'CLim',[0 max(value_array,[],'all')]),ax)
        end
        
        function plot_cluster_shifts(rg, value_array, group_labels_comb)
            arguments
                rg RecordingGroup
                value_array double
                group_labels_comb string
            end
            N_clust = size(value_array,3);
            clust_comp = combvec(1:N_clust,1:N_clust);
            clust_comp = clust_comp(:,clust_comp(1,:) > clust_comp(2,:));
            N_data = sum(sum(value_array(:,:,1)));
            comp_array = zeros(size(value_array));
            
            for i = 1:size(clust_comp,2)
                comp_data = normalize(value_array(:,:,clust_comp(2,i)),'range',[0 1]) - normalize(value_array(:,:,clust_comp(1,i)),'range',[0 1]);
                comp_array(:,:,i) = comp_data./N_data;
                ax(i) = nexttile;
                imagesc(comp_array(:,:,i))
                title(group_labels_comb(clust_comp(2,i)) + " -> " + group_labels_comb(clust_comp(1,i)))
                xticks([])
                yticks([])
            end
            max_change = max(abs(comp_array),[],'all');
            cmap = othercolor('RdBu9',100);
            colormap(cmap)
            arrayfun(@(x) set(x,'CLim',[-max_change max_change]),ax)
        end
        
        function plot_regression_results(rg, regression_var, color_var, color_order)
            arguments
                rg RecordingGroup
                regression_var string
                color_var string = "Mutation"
                color_order double = []
            end
            y_pred = vertcat(rg.Regression.(regression_var).Y_pred);
            y_test = vertcat(rg.Regression.(regression_var).Y_test);
            reg_objects = [rg.Regression.(regression_var).objects];
            [color_idx, color_labels] = rg.returnMetadataArray(reg_objects, color_var);
            
            min_diff = min(diff(unique(y_test)));
            x_vals = y_test/min_diff;
            
            figure('Color','w');
            b = boxchart(x_vals, y_pred,'GroupByColor',color_idx,'MarkerSize',1);
            hold on
            p = plot([0 max(x_vals)*1.2],[0 max(y_test)*1.2],'k--'); p.Color(4) = 0.3;
            xticks(unique(x_vals))
            xticklabels(unique(y_test))
            xlabel(regression_var)
            ylim([0 max(y_test) + 0.5*min(y_test)])
            yticks(unique(y_test))
            ylabel(" Predicted " + regression_var)
            set(gca,'FontSize',8)
            leg = legend(color_labels);
            leg.Box = 'off';
            leg.Location = 'best';
            if isempty(color_order)
               color_order = othercolor('RdBu9',length(b)); 
            end
            arrayfun(@(x) set(b(x),'BoxFaceColor', color_order(x,:)),1:length(b))
        end
        
        function rel_counts = plot_cluster_proportions(rg, separation_var, pooling_vals, clust_method, colors)
            arguments
                rg RecordingGroup
                separation_var string
                pooling_vals cell = {}
                clust_method string = "spectral"
                colors double = []
            end
            [idx, labels] = rg.combineMetadataIndices(rg.DimensionalityReduction.Unit.ObjectGroup,separation_var, pooling_vals);
            cluster_idx = [rg.Clustering.Unit.(clust_method).Index];
            
            cluster_counts_mat = nan(length(labels),length(unique(cluster_idx)));
            for i = 1:length(labels)
                cluster_counts_mat(i,:) = histcounts(cluster_idx(idx == i),'BinLimits',[1 max(cluster_idx)]);
            end
            
            rel_counts = cluster_counts_mat./sum(cluster_counts_mat,2);
%             rel_counts = rel_counts(:,end:-1:1);
            
%             figure('Color','w');
            b = bar(rel_counts,'stacked','FaceColor','flat');
            if ~isempty(colors)
                arrayfun(@(x) set(b(x),'CData',repmat(colors(x,:),[length(labels) 1])),1:length(b))
            end
            xlabel(separation_var)
            xticklabels(labels)
            ylabel("Cluster percentage")
            box off
        end
        
        function plot_group_feature_heatmap(rg, feature_names, separation_var, pooling_vals)
            arguments
                rg RecordingGroup
                feature_names string
                separation_var string
                pooling_vals cell = {}
                
            end
            % TODO
        end

        % function 
        
        function plot_unit_cluster_features(rg, cluster_idx, feature_names, colors)
            % Plots box plots per cluster for the selected features
            arguments
                rg RecordingGroup
                cluster_idx double
                feature_names string
                colors double = []
            end
            input_table = rg.Recordings.getUnitFeatures("all");
            figure('Color','w');
            tiledlayout('flow')
            for f = 1:length(feature_names)
                x = input_table(:,feature_names(f)).Variables;
                nexttile
                b = boxchart(x,'GroupByColor',cluster_idx,'MarkerSize',1);
                title(feature_names(f))
                if ~isempty(colors)
                    arrayfun(@(x) set(b(x), 'BoxFaceColor',colors(x,:),'MarkerColor',colors(x,:)),1:length(b))
                end
                ylim([quantile(x,0.05),quantile(x,0.95)])
            end
        end
        
        function color_mat = plot_unit_cluster_heatmap(rg, cluster_idx, feature_names, cmap)
            % Plots box plots per cluster for the selected features
            arguments
                rg RecordingGroup
                cluster_idx double
                feature_names string = []
                cmap double = othercolor('RdBu9',100);
            end
            N_Clust = length(unique(cluster_idx));
            
            input_table = rg.Recordings.getUnitFeatures("all");
            if isempty(feature_names)
                feature_names = string([input_table.Properties.VariableNames]);
            end
            color_mat = nan(length(feature_names), N_Clust);
            feat_mat = normalize(input_table(:,feature_names).Variables);
            for c = 1:N_Clust
                color_mat(:,c) = mean(feat_mat(cluster_idx == c,:),'omitnan');
            end
            figure('Color','w');
            imagesc(color_mat)
            colormap(cmap)
            xticks(1:N_Clust)
            xticklabels(1:N_Clust)
            xlabel('Cluster ID')
            yticks(1:length(feature_names))
            yticklabels(feature_names)
            colorbar
            set(gca,'CLim',[-max(max(abs(color_mat))) max(max(abs(color_mat)))])
        end
    end
end