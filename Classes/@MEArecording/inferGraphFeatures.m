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
%                  rich_club = rich_club_bd(con_mat,5);
                 nw_feat = ["Density","Assortativity"];%,"RichClub" + (1:5)];
                 nw_val = [density, assortativity];%, rich_club];

             elseif isfield(obj.Connectivity.(alg(a)),"bu") %Binary undirected graphs
                 con_mat = obj.Connectivity.(alg(a)).bu;
                 density = density_und(con_mat);
%                  XY = obj.RecordingInfo.ElectrodeCoordinates([obj.Units.ReferenceElectrode],:);
%                  n = 1000;
%                  tol = 1e-6;
%                  [N,E] = rentian_scaling_2d(con_mat,XY,n,tol);
%                  rm_idx = N == 0 | E == 0;
%                  N(rm_idx) = []; E(rm_idx) = [];
%                  try
%                      [b,~] = robustfit(log10(N),log10(E));
%                      rents_exponent = b(2,1);
%                  catch
rents_exponent = 0;
%                  end
                 clustering_coef = clustering_coef_bu(con_mat);
                 assortativity = assortativity_bin(con_mat,0);
%                  rich_club = rich_club_bu(con_mat,5);
                 nw_feat = ["Density","RentExponent","Assortativity"];%,"RichClub" + (1:5)];
                 nw_val = [density, rents_exponent, assortativity];%, rich_club];

             else
                 disp("No binary connectivity matrix available")
             end

             global_efficiency = efficiency_bin(con_mat);
             local_efficiency = efficiency_bin(con_mat,1);
%              [s,S] = motif3struct_bin(con_mat);
%              s = s./nchoosek(length(con_mat),3); %Normalize
%              S = S./(nchoosek(length(con_mat),3) - nchoosek(length(con_mat) - 1,3)); %Normalize
%              [f,F] = motif3funct_bin(con_mat);
%              f = f./nchoosek(length(con_mat),3); %Normalize
%              F = F./(nchoosek(length(con_mat),3) - nchoosek(length(con_mat) - 1,3)); %Normalize
             eigen_centrality = eigenvector_centrality_und(double(con_mat));
             betweenness = betweenness_bin(con_mat);
             nw_feat = [nw_feat, "GlobalEfficiency"] + "_" + alg(a);%, "StructuralMotif" + (1:13), "FunctionalMotif" + (1:13)] + "_" + alg(a);
             nw_val = [nw_val, global_efficiency];%, s', f'];
             nw_val(isnan(nw_val) | isinf(nw_val)) = 0;
             nw_table = array2table(nw_val,"VariableNames",nw_feat);
             unit_feat = ["ClusteringCoefficient","LocalEfficiency","EigenCentrality","Betweenness"] + "_" + alg(a);%,"UnitStructuralMotif" + (1:13), "UnitFunctionalMotif" + (1:13)];
             unit_val = [clustering_coef, local_efficiency, eigen_centrality, betweenness'];%, S', F'];
             unit_val(isnan(unit_val) | isinf(unit_val)) = 0;
             unit_table = array2table(unit_val,"VariableNames",unit_feat);
             full_nw_table = [full_nw_table nw_table];
             full_unit_table = [full_unit_table unit_table];
          end
          obj.NetworkFeatures.GraphFeatures = full_nw_table;
          for u = 1:size(full_unit_table,1)
              obj.Units(u).GraphFeatures = full_unit_table(u,:);
          end
       end
