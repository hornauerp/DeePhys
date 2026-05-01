function [network_table, unit_table] = computeGraphFeatures(connectivity, graph_params)
% COMPUTEGRAPHFEATURES  Compute graph-theoretic features from connectivity matrices.
%
%   [network_table, unit_table] = computeGraphFeatures(connectivity, graph_params)
%
% INPUTS:
%   connectivity - struct with algorithm subfields (e.g. .CCG, .STTC), each containing
%                  .bd (binary directed) or .bu (binary undirected) connectivity matrices
%   graph_params - struct: .SmallWorldnessIterations (default 100)
%
% OUTPUTS:
%   network_table - (1 x F) table of network-level graph features
%   unit_table    - (N_units x F) table of per-unit graph features

    arguments
        connectivity struct
        graph_params struct = struct()
    end

    if ~isfield(graph_params, 'SmallWorldnessIterations')
        graph_params.SmallWorldnessIterations = 100;
    end

    alg = string(fieldnames(connectivity));
    full_nw_table = table();
    full_unit_table = table();

    for a = 1:length(alg)
        if isfield(connectivity.(alg(a)), 'bd')
            con_mat = connectivity.(alg(a)).bd;
            density = density_dir(con_mat);
            clustering_coef = clustering_coef_bd(con_mat);
            assortativity = assortativity_bin(con_mat, 1);
            rich_club = rich_club_bd(con_mat, min(5, max(sum(con_mat ~= 0))));
            nw_feat = ["Density","Assortativity","RichClub" + (1:length(rich_club))];
            nw_val = [density, assortativity, rich_club];
            is_directed = true;

        elseif isfield(connectivity.(alg(a)), 'bu')
            con_mat = connectivity.(alg(a)).bu;
            density = density_und(con_mat);
            clustering_coef = clustering_coef_bu(con_mat);
            assortativity = assortativity_bin(con_mat, 0);
            rich_club = rich_club_bu(con_mat, min(5, max(sum(con_mat))));
            nw_feat = ["Density","Assortativity","RichClub" + (1:length(rich_club))];
            nw_val = [density, assortativity, rich_club];
            is_directed = false;

        else
            warning('computeGraphFeatures:noMatrix', ...
                'Algorithm "%s" has no binary connectivity matrix — skipping.', alg(a));
            continue
        end

        global_efficiency = efficiency_bin(con_mat);
        local_efficiency = efficiency_bin(con_mat, 1);

        % community_louvain (BCT) uses randperm internally, making each run
        % non-deterministic. Running multiple iterations and taking max(Q)
        % reduces variance. n_louvain=5 balances stability vs. compute cost.
        n_louvain = 5;
        Q_vals = nan(1, n_louvain);
        for lr = 1:n_louvain
            try
                [~, Q_vals(lr)] = community_louvain(con_mat);
            catch
            end
        end
        Q_modularity = max(Q_vals, [], 'omitnan');
        if isempty(Q_modularity) || all(isnan(Q_vals))
            Q_modularity = NaN;
        end

        n_rand = graph_params.SmallWorldnessIterations;
        n_nodes = size(con_mat, 1);
        if density > 0 && n_nodes >= 4
            C_real = mean(clustering_coef);
            L_real = charpath(distance_bin(con_mat));
            C_rand_vals = zeros(1, n_rand);
            L_rand_vals = zeros(1, n_rand);
            for ir = 1:n_rand
                if is_directed
                    rand_mat = randmio_dir(double(con_mat ~= 0), 5);
                    C_rand_vals(ir) = mean(clustering_coef_bd(rand_mat));
                else
                    rand_mat = randmio_und(double(con_mat ~= 0), 5);
                    C_rand_vals(ir) = mean(clustering_coef_bu(rand_mat));
                end
                L_rand_vals(ir) = charpath(distance_bin(rand_mat));
            end
            C_rand = mean(C_rand_vals);
            L_rand = mean(L_rand_vals);
            if C_rand > 0 && L_rand > 0
                small_worldness = (C_real / C_rand) / (L_real / L_rand);
            else
                small_worldness = NaN;
            end
        else
            small_worldness = NaN;
        end

        if is_directed
            eigen_centrality = eigenvector_centrality_dir(double(con_mat ~= 0));
        else
            eigen_centrality = eigenvector_centrality_und(double(con_mat ~= 0));
        end
        betweenness = betweenness_bin(con_mat);

        nw_feat = [nw_feat, "GlobalEfficiency", "Modularity", "SmallWorldness"] + "_" + alg(a);
        nw_val = [nw_val, global_efficiency, Q_modularity, small_worldness];
        nw_val(isinf(nw_val)) = NaN;
        nw_table = array2table(nw_val, "VariableNames", nw_feat);

        unit_feat = ["ClusteringCoefficient","LocalEfficiency","EigenCentrality","Betweenness"] + "_" + alg(a);
        unit_val = [clustering_coef, local_efficiency, eigen_centrality, betweenness'];
        unit_val(isinf(unit_val)) = NaN;
        unit_table_alg = array2table(unit_val, "VariableNames", unit_feat);

        full_nw_table = [full_nw_table nw_table]; %#ok<AGROW>
        full_unit_table = [full_unit_table unit_table_alg]; %#ok<AGROW>
    end

    network_table = full_nw_table;
    unit_table = full_unit_table;
end
