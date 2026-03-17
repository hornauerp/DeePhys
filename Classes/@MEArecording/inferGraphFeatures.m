       function inferGraphFeatures(obj, alg)
          arguments
             obj MEArecording
             alg (1,:) string = [] %Connectivity inference algorithms to compute graphs for
          end

          if isempty(alg) %Compute features for all available inferred graphs
              alg = string(fields(obj.Connectivity));
              assert(~isempty(alg),"No connectivity inference results found")
          end

          full_nw_table = table();
          full_unit_table = table();
          for a = 1:length(alg)
             assert(isfield(obj.Connectivity,alg(a)),"No connectivity inference results found")
             if isfield(obj.Connectivity.(alg(a)),"bd") %Binary directed graphs
                 con_mat = obj.Connectivity.(alg(a)).bd;
                 density = density_dir(con_mat);
                 clustering_coef = clustering_coef_bd(con_mat);
                 assortativity = assortativity_bin(con_mat,1);
                 rich_club = rich_club_bd(con_mat, min(5, max(sum(con_mat ~= 0))));
                 nw_feat = ["Density","Assortativity","RichClub" + (1:length(rich_club))];
                 nw_val = [density, assortativity, rich_club];

             elseif isfield(obj.Connectivity.(alg(a)),"bu") %Binary undirected graphs
                 con_mat = obj.Connectivity.(alg(a)).bu;
                 density = density_und(con_mat);
                 clustering_coef = clustering_coef_bu(con_mat);
                 assortativity = assortativity_bin(con_mat,0);
                 rich_club = rich_club_bu(con_mat, min(5, max(sum(con_mat))));
                 nw_feat = ["Density","Assortativity","RichClub" + (1:length(rich_club))];
                 nw_val = [density, assortativity, rich_club];

             else
                 warning('MEArecording:inferGraphFeatures', ...
                     'Algorithm "%s" has no binary connectivity matrix — skipping.', alg(a));
                 continue
             end

             global_efficiency = efficiency_bin(con_mat);
             local_efficiency = efficiency_bin(con_mat,1);

             % Modularity (Louvain community detection)
             try
                 [~, Q_modularity] = community_louvain(con_mat);
             catch
                 Q_modularity = NaN;
             end

             % Small-worldness: sigma = (C/C_rand) / (L/L_rand)
             % Using 20 random rewirings for comparison
             n_rand = 20;
             n_nodes = size(con_mat, 1);
             is_directed = isfield(obj.Connectivity.(alg(a)), "bd");
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

             eigen_centrality = eigenvector_centrality_und(double(con_mat ~= 0));
             betweenness = betweenness_bin(con_mat);
             nw_feat = [nw_feat, "GlobalEfficiency", "Modularity", "SmallWorldness"] + "_" + alg(a);
             nw_val = [nw_val, global_efficiency, Q_modularity, small_worldness];
             nw_val(isinf(nw_val)) = NaN;
             nw_table = array2table(nw_val,"VariableNames",nw_feat);
             unit_feat = ["ClusteringCoefficient","LocalEfficiency","EigenCentrality","Betweenness"] + "_" + alg(a);%,"UnitStructuralMotif" + (1:13), "UnitFunctionalMotif" + (1:13)];
             unit_val = [clustering_coef, local_efficiency, eigen_centrality, betweenness'];%, S', F'];
             unit_val(isinf(unit_val)) = NaN;
             unit_table = array2table(unit_val,"VariableNames",unit_feat);
             full_nw_table = [full_nw_table nw_table];
             full_unit_table = [full_unit_table unit_table];
          end
          obj.NetworkFeatures.GraphFeatures = full_nw_table;
          for u = 1:size(full_unit_table,1)
              obj.Units(u).GraphFeatures = full_unit_table(u,:);
          end
       end
