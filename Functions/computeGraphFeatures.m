function [network_table, unit_table] = computeGraphFeatures(connectivity, graph_params)
% COMPUTEGRAPHFEATURES  Compute graph-theoretic features from connectivity matrices.
%
%   [network_table, unit_table] = computeGraphFeatures(connectivity, graph_params)
%
% INPUTS:
%   connectivity - struct with algorithm subfields (e.g. .CCG, .STTC, .DDC), each containing
%                  .bd (binary directed) or .bu (binary undirected) connectivity matrices,
%                  and optionally link weights (.wd_sig / .wd for directed, .wu for
%                  undirected graphs)
%   graph_params - struct: .SmallWorldnessIterations (default 100)
%                          .Types: "binary" (default), "weighted" or both. Weighted
%                          features use |weights| of the significant links, normalised to
%                          (0, 1], and are named <Feature>_<Algorithm>_w. Algorithms
%                          without weights (CCG) only give binary features.
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
    types = "binary";
    if isfield(graph_params, 'Types') && ~isempty(graph_params.Types)
        types = lower(string(graph_params.Types));
    end
    do_binary   = any(types == "binary");
    do_weighted = any(types == "weighted");

    alg = string(fieldnames(connectivity));
    full_nw_table = table();
    full_unit_table = table();

    for a = 1:length(alg)
        if do_weighted
            [nw_w, unit_w] = weightedGraphFeatures(connectivity.(alg(a)), alg(a));
            if ~isempty(nw_w)
                full_nw_table   = [full_nw_table nw_w]; %#ok<AGROW>
                full_unit_table = [full_unit_table unit_w]; %#ok<AGROW>
            end
        end
        if ~do_binary
            continue
        end
        if isfield(connectivity.(alg(a)), 'bd')
            con_mat = double(connectivity.(alg(a)).bd ~= 0);   % BCT binary measures need 0/1; CCG codes inhibitory links as -1
            density = density_dir(con_mat);
            clustering_coef = clustering_coef_bd(con_mat);
            assortativity = assortativity_bin(con_mat, 1);
            rich_club = richClubLevels(@rich_club_bd, con_mat, max(sum(con_mat ~= 0)));
            nw_feat = ["Density","Assortativity","RichClub" + (1:length(rich_club))];
            nw_val = [density, assortativity, rich_club];
            is_directed = true;

        elseif isfield(connectivity.(alg(a)), 'bu')
            con_mat = double(connectivity.(alg(a)).bu ~= 0);
            density = density_und(con_mat);
            clustering_coef = clustering_coef_bu(con_mat);
            assortativity = assortativity_bin(con_mat, 0);
            rich_club = richClubLevels(@rich_club_bu, con_mat, max(sum(con_mat)));
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
        % reduces variance. n_louvain=20 gives stable estimates with modest cost.
        n_louvain = 20;
        Q_vals = nan(1, n_louvain);
        for lr = 1:n_louvain
            try
                [~, Q_vals(lr)] = community_louvain(con_mat);
            catch ME
                warning('computeGraphFeatures:louvainFailed', ...
                    'community_louvain failed on replicate %d: %s', lr, ME.message);
            end
        end
        Q_modularity = max(Q_vals, [], 'omitnan');
        if isempty(Q_modularity) || all(isnan(Q_vals))
            Q_modularity = NaN;
        end

        n_rand = graph_params.SmallWorldnessIterations;
        n_nodes = size(con_mat, 1);
        % Skip small-worldness for sparse networks. Path lengths are averaged over reachable
        % pairs only, so one isolated node does not make L infinite (and SW NaN or 0).
        % randmio_dir / randmio_und loop forever unless two edges have four distinct
        % end nodes (e.g. a graph whose only links are one reciprocal pair)
        if density > 0 && density >= 0.05 && n_nodes >= 4 && hasDisjointEdgePair(con_mat)
            C_real = mean(clustering_coef);
            L_real = charpath(distance_bin(con_mat), 0, 0);   % finite paths only (BCT default counts Inf)
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
                L_rand_vals(ir) = charpath(distance_bin(rand_mat), 0, 0);
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

function [nw_table, unit_table] = weightedGraphFeatures(c, name)
% Weighted graph of the significant links: |link weights|, normalised to (0, 1].
% Directed: DDC (.wd_sig, or .wd masked by .bd). Undirected: STTC (.wu masked by .bu;
% negative tiling coefficients are no links). Returns empty tables without weights.
    nw_table = table(); unit_table = table();
    if isfield(c, 'wd_sig')
        W = abs(c.wd_sig); directed = true;
    elseif isfield(c, 'wd') && isfield(c, 'bd')
        W = abs(c.wd) .* (c.bd ~= 0); directed = true;
    elseif isfield(c, 'wu') && isfield(c, 'bu')
        W = max(c.wu, 0) .* (c.bu ~= 0); directed = false;
    else
        return
    end
    nw_names   = ["MeanStrength", "Assortativity", "GlobalEfficiency", "Modularity"] + "_" + name + "_w";
    unit_names = ["Strength", "ClusteringCoefficient", "LocalEfficiency", "EigenCentrality", ...
                  "Betweenness"] + "_" + name + "_w";
    n = size(W, 1);
    W(~isfinite(W)) = 0;
    W(1:n+1:end) = 0;
    if ~any(W(:))   % no significant links: features undefined
        nw_table   = array2table(nan(1, numel(nw_names)), 'VariableNames', nw_names);
        unit_table = array2table(nan(n, numel(unit_names)), 'VariableNames', unit_names);
        return
    end
    W = W / max(W(:));

    if directed
        [~, ~, strength] = strengths_dir(W);       % in + out strength
        clustering = clustering_coef_wd(W);
        assort     = assortativity_wei(W, 1);      % out-strength / in-strength
        eig_c      = eigenvector_centrality_dir(W);
    else
        strength   = strengths_und(W);
        clustering = clustering_coef_wu(W);
        assort     = assortativity_wei(W, 0);
        eig_c      = eigenvector_centrality_und(W);
    end
    glob_eff = efficiency_wei(W);
    loc_eff  = nan(n, 1);
    try
        loc_eff = efficiency_wei(W, 2);            % modified local efficiency (Wang et al. 2017)
    catch
    end
    betw = betweenness_wei(weight_conversion(W, 'lengths'));   % paths along strong links
    Q = NaN;
    for r = 1:20
        try
            [~, q] = community_louvain(W);
            Q = max(Q, q);
        catch
        end
    end
    nw_val = [mean(strength), assort, glob_eff, Q];
    nw_val(isinf(nw_val)) = NaN;
    unit_val = [strength(:), clustering(:), loc_eff(:), eig_c(:), betw(:)];
    unit_val(isinf(unit_val)) = NaN;
    nw_table   = array2table(nw_val, 'VariableNames', nw_names);
    unit_table = array2table(unit_val, 'VariableNames', unit_names);
end

function rc = richClubLevels(rc_fn, C, max_degree)
% Rich-club coefficients at degree levels 1..5, NaN above the maximum degree: a fixed set of
% columns, so recordings with different maximum degrees stack into one table.
    n_levels = 5;
    rc = nan(1, n_levels);
    k = min(n_levels, max_degree);
    if k >= 1
        v = rc_fn(C, k);
        rc(1:numel(v)) = v;
    end
end

function tf = hasDisjointEdgePair(C)
% True if two links of C share no end node (needed by the BCT rewiring randomisers).
    [i, j] = find(triu(C | C', 1));   % undirected view: a reciprocal pair is one link
    tf = false;
    for e = 1:numel(i)
        if any(i ~= i(e) & i ~= j(e) & j ~= i(e) & j ~= j(e))
            tf = true;
            return
        end
    end
end
