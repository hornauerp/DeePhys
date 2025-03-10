addpath(genpath("/home/phornauer/Git/DeePhys"))
set(0,'defaultFigureRenderer','painters')

%% LOAD DATA
save_root = "/home/phornauer/Git/DeePhys/Data/Supplement/";
load(fullfile(save_root,'figure_S3.mat'))

%% SET PARAMETERS
plot_centroid = false;
nodeSz = 10;
mapSz = 300;
sigma = mapSz/80;
composition_cmap = othercolor('YlGn5',2);
mutation_cmap = othercolor('Spectral9',2);
fontsz = 6;

fs3 = figure('Color','w','Position',[1200 200 595 595],'Units','points');
tl = tiledlayout(5,4,'TileSpacing','compact','Padding','tight');

nt(1) = nexttile(tl,1);
RecordingGroup.plot_cluster_outlines(composition_reduction, composition_idx, nt(1), plot_centroid, nodeSz, mapSz, sigma, composition_cmap);
leg(1) = legend(nt(1),["Homogeneous", "Heterogeneous"],'Box','off');
leg(1).Position(1:2) = nt(1).Position(1:2);
leg(1).Position(2) = leg(1).Position(2) - 0.1;
leg(1).ItemTokenSize = [4 4];

nt(2) = nexttile(tl,2);
RecordingGroup.plot_cluster_outlines(mutation_reduction, mutation_idx, nt(2), plot_centroid, nodeSz, mapSz, sigma, mutation_cmap);
leg(2) = legend(nt(2),["Control", "Mutant"],'Box','off');
leg(2).Position(1:2) = nt(2).Position(1:2);
leg(2).Position(2) = leg(2).Position(2) - 0.1;
leg(2).ItemTokenSize = [4 4];
setAllFontSizes(fs3, fontsz)

%%
exportgraphics(fs3,'/home/phornauer/Git/DeePhys/Plots/FigureS3.tif','Resolution',300)