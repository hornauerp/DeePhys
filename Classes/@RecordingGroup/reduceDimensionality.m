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

            t_start = tic;

            if isempty(grouping_var) %Get group indices for stratification
                object_group = [rg.Recordings]; %Stratify on recordings level also for units, to avoid bias
            else
                [input_table, object_group] = aggregateCultureFeatureTables(rg, level, grouping_var, grouping_values, tolerance, network_features, unit_features, normalization, useClustered);
            end

            if isempty(grouping_var)
                if level == "Unit"
                    input_table = MEArecording.getUnitFeatures(object_group,unit_features);
                    object_group = [rg.Units];
                    n_neighbors = 100;
                elseif level == "Recording"
                    input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
                    n_neighbors = 100;
                else
                    error('Unknown level')
                end
            else
                if level == "Unit"
                    object_group = cellfun(@(x) [x(1).Units],object_group,'un',0);
                    object_group = [object_group{:}];
                    n_neighbors = 100;
                else
                    n_neighbors = 10;
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

                    color_file = fullfile(rg.Recordings(1).getParentPath(),'umap','colorsByName.properties');
                    [reduction, umap, clusterIdentifiers, extras] = run_umap(norm_data,'n_components',n_dims,'n_neighbors',n_neighbors,'min_dist',1,'cluster_detail','adaptive','spread',5,'sgd_tasks',20,...
                        'verbose','none','color_file',color_file); %Torsten
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

            % Build result object
            dr_result = DimReductionResult( ...
                'Reduction', reduction, ...
                'ObjectGroup', object_group, ...
                'Method', method, ...
                'UnitFeatures', unit_features, ...
                'NetworkFeatures', network_features, ...
                'Parameters', struct('n_dims', n_dims, 'normalization_var', normalization_var, ...
                    'grouping_var', grouping_var, 'grouping_values', grouping_values));

            if method == "UMAP"
                dr_result.Graph = rg.DimensionalityReduction.(level).(method).Graph;
            end

            % Backward compat: still store on rg
            rg.DimensionalityReduction.(level).(method).Reduction = reduction;
            rg.DimensionalityReduction.(level).(method).GroupingVariable = grouping_var;
            rg.DimensionalityReduction.(level).(method).UnitFeatures = unit_features;
            rg.DimensionalityReduction.(level).(method).NetworkFeatures = network_features;
            rg.DimensionalityReduction.(level).ObjectGroup = object_group;

            % Log the analysis
            elapsed = toc(t_start);
            AnalysisLog.instance().add('reduceDimensionality', ...
                struct('level', level, 'method', method, 'n_dims', n_dims), ...
                sprintf('%d objects, %s', size(reduction, 1), method), elapsed);
end
