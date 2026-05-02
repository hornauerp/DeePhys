function [network_table, unit_table] = computeCellTypeGraphFeatures(connectivity, cell_type_labels, graph_params)
% COMPUTECELLTYPEGRAPHFEATURES  Graph features on E/I sub-networks.
%
%   [network_table, unit_table] = computeCellTypeGraphFeatures(connectivity, cell_type_labels, graph_params)
%
% For each connectivity algorithm, extracts excitatory-only (EE) and
% inhibitory-only (II) sub-matrices and delegates to computeGraphFeatures.
% Cross-type connection density is computed directly (rectangular sub-matrices
% are not valid inputs for BCT graph functions).
%
% INPUTS:
%   connectivity     - struct with algorithm subfields (.CCG, .STTC, .DDC),
%                      each containing .bd (binary directed) or .bu (binary undirected)
%   cell_type_labels - (1 x N) double: 1=excitatory, 2=inhibitory, NaN=unclassified
%   graph_params     - struct passed through to computeGraphFeatures
%
% OUTPUTS:
%   network_table - (1 x F) table of network-level cell-type graph features
%   unit_table    - (N x F) table of per-unit cell-type features

    arguments
        connectivity     struct
        cell_type_labels (1,:) double
        graph_params     struct = struct()
    end

    n_units = numel(cell_type_labels);
    exc_idx = find(cell_type_labels == 1);
    inh_idx = find(cell_type_labels == 2);

    alg = string(fieldnames(connectivity));
    full_nw_table = table();
    full_unit_table = table();

    for a = 1:numel(alg)
        alg_data = connectivity.(alg(a));

        % Determine matrix type
        if isfield(alg_data, 'bd')
            con_mat = alg_data.bd;
            mat_field = "bd";
        elseif isfield(alg_data, 'bu')
            con_mat = alg_data.bu;
            mat_field = "bu";
        else
            continue
        end

        % --- Square sub-networks (EE, II) via computeGraphFeatures ---
        sub_defs = struct('name', {"EE", "II"}, 'idx', {exc_idx, inh_idx});
        for s = 1:numel(sub_defs)
            idx = sub_defs(s).idx;
            suffix = "_" + sub_defs(s).name;

            if numel(idx) >= 3
                sub_mat = con_mat(idx, idx);
                sub_conn = struct();
                sub_conn.(alg(a)).(mat_field) = sub_mat;
                [nw_sub, ~] = computeGraphFeatures(sub_conn, graph_params);
                nw_sub.Properties.VariableNames = nw_sub.Properties.VariableNames + suffix;
            else
                % Degenerate: fewer than 3 nodes — produce NaN placeholders
                nw_sub = nanPlaceholderNetwork(alg(a), mat_field, suffix, graph_params);
            end
            full_nw_table = [full_nw_table, nw_sub]; %#ok<AGROW>
        end

        % --- Cross-type density ---
        is_directed = mat_field == "bd";
        if ~isempty(exc_idx) && ~isempty(inh_idx)
            if is_directed
                ei_mat = con_mat(exc_idx, inh_idx);
                ie_mat = con_mat(inh_idx, exc_idx);
                ei_density = nnz(ei_mat) / numel(ei_mat);
                ie_density = nnz(ie_mat) / numel(ie_mat);
                ct_names = ["CrossTypeDensity_" + alg(a) + "_EI", ...
                            "CrossTypeDensity_" + alg(a) + "_IE"];
                ct_vals = [ei_density, ie_density];
            else
                cross_mat = con_mat(exc_idx, inh_idx);
                cross_density = nnz(cross_mat) / numel(cross_mat);
                ct_names = "CrossTypeDensity_" + alg(a);
                ct_vals = cross_density;
            end
        else
            if is_directed
                ct_names = ["CrossTypeDensity_" + alg(a) + "_EI", ...
                            "CrossTypeDensity_" + alg(a) + "_IE"];
                ct_vals = [NaN, NaN];
            else
                ct_names = "CrossTypeDensity_" + alg(a);
                ct_vals = NaN;
            end
        end
        ct_table = array2table(ct_vals, 'VariableNames', ct_names);
        full_nw_table = [full_nw_table, ct_table]; %#ok<AGROW>

        % --- Per-unit within-type and cross-type connection counts ---
        within_conn = NaN(n_units, 1);
        cross_conn  = NaN(n_units, 1);

        for u = 1:n_units
            lbl = cell_type_labels(u);
            if isnan(lbl), continue; end
            if lbl == 1
                same_idx = exc_idx;
                other_idx = inh_idx;
            else
                same_idx = inh_idx;
                other_idx = exc_idx;
            end

            if is_directed
                out_row = con_mat(u, :);
                in_col  = con_mat(:, u)';
                within_conn(u) = sum(abs(out_row(same_idx)))  + sum(abs(in_col(same_idx)));
                cross_conn(u)  = sum(abs(out_row(other_idx))) + sum(abs(in_col(other_idx)));
            else
                row = con_mat(u, :);
                within_conn(u) = sum(abs(row(same_idx)));
                cross_conn(u)  = sum(abs(row(other_idx)));
            end
        end

        unit_alg = table(within_conn, cross_conn, ...
            'VariableNames', ["WithinTypeConnections_" + alg(a), ...
                              "CrossTypeConnections_" + alg(a)]);
        full_unit_table = [full_unit_table, unit_alg]; %#ok<AGROW>
    end

    network_table = full_nw_table;
    unit_table    = full_unit_table;
end


function nw_table = nanPlaceholderNetwork(alg_name, mat_field, suffix, graph_params)
% Build a one-row NaN table matching computeGraphFeatures output column names.
    if ~isfield(graph_params, 'SmallWorldnessIterations')
        graph_params.SmallWorldnessIterations = 100;
    end

    % Match the column naming from computeGraphFeatures
    n_rich_club = 5;
    nw_feat = ["Density", "Assortativity", "RichClub" + (1:n_rich_club), ...
               "GlobalEfficiency", "Modularity", "SmallWorldness"] + "_" + alg_name + suffix;
    nw_val = NaN(1, numel(nw_feat));
    nw_table = array2table(nw_val, 'VariableNames', nw_feat);
end
