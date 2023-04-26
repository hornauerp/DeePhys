addpath(genpath("/home/phornauer/Git/DeePhys"))
set(0,'defaultFigureRenderer','painters')

%% LOAD DATA
save_root = "/home/phornauer/Git/DeePhys/Data/Figure6";
load(fullfile(save_root,'qp_results.mat'))

%% PREPARE DATA
outlier_idx = qp_results.response.Unit.umap(:,1) > 10;
qp_results.response.Unit.cluster_idx(outlier_idx) = [];
qp_results.response.Unit.umap(outlier_idx,:) = [];
qp_results.response.Unit.group_idx(outlier_idx) = [];
qp_results.baseline.Unit.cluster_idx(outlier_idx) = [];
qp_results.baseline.Unit.umap(outlier_idx,:) = [];

jitter = linspace(-0.1,0.1,size(qp_results.response.mean,3));
x_timeline = 1:size(qp_results.response.mean,2);
x_mat = x_timeline([ones(size(jitter))],:) - jitter';
xx_timeline = linspace(min(x_timeline),max(x_timeline),100);
xx_mat = xx_timeline([ones(size(jitter))],:) - jitter';
y{1} = permute(qp_results.response.mean(26,:,:),[3,2,1]);
y{2} = permute(qp_results.response.mean(27,:,:),[3,2,1]);
yy_timeline = cellfun(@(i) pchip(x_timeline,i,xx_timeline),y,'un',0);
sd{1} = permute(qp_results.response.sd(26,:,:),[3,2,1]);
sd{2} = permute(qp_results.response.sd(27,:,:),[3,2,1]);

%% SET PARAMETERS
greens = othercolor('Greens9',100); %greens = greens(30:end,:);
plot_centroid = false;
sc_nodeSz = 3;
nw_nodeSz = 10;
mapSz = 300;
sigma = mapSz/80;
qp_cmap = othercolor('Mdarkrainbow',max(qp_results.response.Unit.cluster_idx));
bl_cmap = othercolor('Mdarkrainbow',max(qp_results.baseline.Unit.cluster_idx));
% act_cmap = othercolor('Mdarkrainbow',max(act.cluster_idx));
% comb_cmap = othercolor('Mdarkrainbow',max(comb.cluster_idx));
fontsz = 6;

%% NEW VERSION
f6 = figure('Color','w','Position',[1200 200 595 595],'Units','points');
tl = tiledlayout(5,4,'TileSpacing','compact','Padding','tight');

nt(1) = nexttile(tl,1);
RecordingGroup.plot_cluster_outlines(qp_results.response.Unit.umap, qp_results.response.Unit.cluster_idx, nt(1), plot_centroid, sc_nodeSz, mapSz, sigma, qp_cmap)
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(1).Children(1:end-1))
set(nt(1),'XLim',[50 250], 'YLim', [50 250])
title('Single-Cell QP Response','FontWeight','normal')

nt(2) = nexttile(tl,2);
RecordingGroup.plot_cluster_outlines(qp_results.baseline.Unit.umap, ones(1,length(qp_results.baseline.Unit.umap)), nt(2), plot_centroid, sc_nodeSz, mapSz, sigma, [0.5 0.5 0.5])
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(2).Children(1:end-1))
set(nt(2),'XLim',[50 250], 'YLim', [50 250])
title('Single-Cell Baseline','FontWeight','normal')

nt(3) = nexttile(tl,3);
RecordingGroup.plot_cluster_outlines(qp_results.baseline.Unit.umap, qp_results.response.Unit.cluster_idx, nt(3), plot_centroid, sc_nodeSz, mapSz, sigma, qp_cmap)
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(3).Children(1:end-1))
ch = get(nt(3),'Children');
set(nt(3),'XLim',[50 250], 'YLim', [50 250], 'Children',ch([2,1,3,4,5]))
title('QP Response Mapping','FontWeight','normal')

nt(4) = nexttile(tl,4);
RecordingGroup.plot_cluster_outlines(qp_results.response.Network.clustered.umap, qp_results.response.Network.clustered.cluster_idx, nt(4), plot_centroid, nw_nodeSz, mapSz, sigma, qp_results.response.cmap([4,3,2,1],:))
set(nt(4),'XLim',[40 270], 'YLim', [40 270])
leg = legend(nt(4),["SNCA + LNA","SNCA","WT + LNA","WT"],'Location','north','NumColumns',2,'Box','off');
leg.ItemTokenSize = [4 4]; leg.Position(2) = 0.96;
% title(nt(4), 'QP Response Clustered','FontWeight','normal')

setAllFontSizes(f6, fontsz)

%%
exportgraphics(f6,'/home/phornauer/Git/DeePhys/Plots/Figure6.tif','Resolution',300)

%% OLD VERSION
f6 = figure('Color','w','Position',[1200 200 595 595],'Units','points');
tl = tiledlayout(5,4,'TileSpacing','compact','Padding','tight');

% UMAP WF clustering
nt(1) = nexttile(tl,1);
RecordingGroup.plot_cluster_outlines(qp_results.response.Unit.umap, qp_results.response.Unit.cluster_idx, nt(1), plot_centroid, sc_nodeSz, mapSz, sigma, qp_cmap)
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(1).Children(1:end-1))
set(nt(1),'XLim',[50 250], 'YLim', [50 250])
title('Single-Cell QP Response','FontWeight','normal')

nt(2) = nexttile(tl,2);
bar([5,4,2,1],qp_results.response.cluster_proportions,'stacked')
set(nt(2),'Box','off','XTickLabel',["       WT\newline      -","\newline LNA","       SNCA\newline       -","\newline LNA"],'ColorOrder',qp_cmap,'XTickLabelRotation',0)
ylabel('Cluster Proportions'); %title(nt(2),'Activity Classification','FontWeight','normal')

nt(3) = nexttile(tl,3);
plot(nt(3),xx_mat',yy_timeline{1}'); hold on
errorbar(x_mat',y{1}',sd{1}',...
    'LineWidth',1,'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',2,...
    'HandleVisibility','off');
set(nt(3),'ColorOrder',[qp_cmap;qp_cmap],'Box','off','XTickLabel',[0 10 20])
xlabel(['Quinpirole [' char(181) 'M]']); ylabel('FAR')

nt(4) = nexttile(tl,4);
plot(nt(4),xx_mat',yy_timeline{2}'); hold on
errorbar(x_mat',y{2}',sd{2}',...
    'LineWidth',1,'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',2,...
    'HandleVisibility','off');
set(nt(4),'ColorOrder',[qp_cmap; qp_cmap],'Box','off','XTickLabel',[0 10 20])
xlabel(['Quinpirole [' char(181) 'M]']); ylabel('FAD')

nt(5) = nexttile(tl,5);
RecordingGroup.plot_cluster_outlines(qp_results.baseline.Unit.umap, qp_results.baseline.Unit.cluster_idx, nt(5), plot_centroid, sc_nodeSz, mapSz, sigma, bl_cmap)
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(5).Children(1:end-1))
set(nt(5),'XLim',[50 250], 'YLim', [50 250])
title('Single-Cell Baseline','FontWeight','normal')

nt(6) = nexttile(tl,6);
RecordingGroup.plot_cluster_outlines(qp_results.baseline.Unit.umap, qp_results.response.Unit.cluster_idx, nt(6), plot_centroid, sc_nodeSz, mapSz, sigma, qp_cmap)
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(6).Children(1:end-1))
ch = get(nt(6),'Children');
set(nt(6),'XLim',[50 250], 'YLim', [50 250], 'Children',ch([2,1,3,4,5]))
title('QP Response Mapping','FontWeight','normal')

nt(7) = nexttile(tl,7);
RecordingGroup.plot_cluster_outlines(qp_results.response.Network.unclustered.umap, qp_results.response.Network.unclustered.cluster_idx, nt(7), plot_centroid, nw_nodeSz, mapSz, sigma, qp_results.response.cmap([4,3,2,1],:))
ch = get(nt(7),'Children'); 
set(nt(7),'XLim',[40 260], 'YLim', [40 260],'Children',ch([4,3,2,1,5]))
title(nt(7), 'QP Response Unclustered','FontWeight','normal')

leg(1) = legend(nt(7),["WT","WT +\newlineLNA", "SNCA", "SNCA +\newlineLNA"],'Location','eastoutside','NumColumns',1,'Box','off');
leg(1).ItemTokenSize = [4 4]; leg(1).Position= leg(1).Position+0.0001;

nt(8) = nexttile(tl,8);
RecordingGroup.plot_cluster_outlines(qp_results.response.Network.clustered.umap, qp_results.response.Network.clustered.cluster_idx, nt(8), plot_centroid, nw_nodeSz, mapSz, sigma, qp_results.response.cmap([4,3,2,1],:))
set(nt(8),'XLim',[40 260], 'YLim', [40 260])
title(nt(8), 'QP Response Clustered','FontWeight','normal')

setAllFontSizes(f6, fontsz)
