function assignUnitClusterIdx(rg, method, calc_feat, concatenated)
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
