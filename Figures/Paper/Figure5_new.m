addpath(genpath("/home/phornauer/Git/DeePhys"))
set(0,'defaultFigureRenderer','painters')

%% LOAD DATA
save_root = "/home/phornauer/Git/DeePhys/Data/Figure5";
load(fullfile(save_root,'activity_prediction.mat'))
load(fullfile(save_root,'waveform_prediction.mat'))
load(fullfile(save_root,'combined_prediction.mat'))
load(fullfile(save_root,'cluster_scatter.mat'))
load(fullfile(save_root,'scatter_comb.mat'))
load(fullfile(save_root,'feature_group_heatmap_pooled_no_lna.mat'))
% load(fullfile(save_root,'single_cell_predictions.mat'))
load(fullfile(save_root,'cluster_proportions.mat'))

%% PREPARE DATA
x_wf = (0:length(wf.pI.mean))/20; % Adjust for missing peak (0 variance)
xx_wf = 0:0.001:length(x_wf); 
wf_pI = [wf.pI.mean(1:19) 0 wf.pI.mean(20:60)];
yy_wf_pI = pchip(x_wf,wf_pI,xx_wf);
yy_wf_pI(yy_wf_pI<0) = 0;
mean_wf_matrix = [wf.mean_wfs{:}]';
interp_wfs = pchip(x_wf,mean_wf_matrix,xx_wf);
interp_mean_wf = mean(interp_wfs);

mean_wf_comb_matrix = [comb.mean_wfs{:}]';
interp_comb_wfs = pchip(x_wf,mean_wf_comb_matrix,xx_wf);
interp_comb_mean_wf = mean(interp_comb_wfs);

norm_cm_wf = wf.confusion_mat./sum(wf.confusion_mat,2);
norm_cm_act = act.confusion_mat./sum(act.confusion_mat,2);
norm_cm_comb = comb.confusion_mat./sum(comb.confusion_mat,2);

norm_wf_pred = wf_table.Variables./sum(wf_table.Variables,2);
norm_act_pred = act_table.Variables./sum(act_table.Variables,2);
norm_comb_pred = comb_table.Variables./sum(comb_table.Variables,2);



%% SET PARAMETERS
greens = othercolor('Greens9',100); %greens = greens(30:end,:);
plot_centroid = false;
nodeSz = 3;
mapSz = 300;
sigma = mapSz/80;
wf_cmap = othercolor('Mdarkrainbow',max(wf.cluster_idx));
act_cmap = othercolor('Mdarkrainbow',max(act.cluster_idx));
comb_cmap = othercolor('Mdarkrainbow',max(comb.cluster_idx));
fontsz = 6.5;
cb_offset = 0.005;

act_features = ["sSFD","sSFR","sAMI","MIS","sRM"];
clust_prop_x = [1:2,4:5,8:9,11:12];
clust_prop_x_labels = ["     -\newline     WT", "+","      -\newline     SNCA", "+", "     -\newline     WT", "+", "      -\newline     A53T", "+"];

%%
f5 = figure('Color','w','Position',[1200 200 595 595],'Units','points');
tl = tiledlayout(5,12,'TileSpacing','tight','Padding','tight');

% UMAP WF clustering
nt(1) = nexttile(tl,1,[1 3]);
RecordingGroup.plot_cluster_outlines(wf.umap, wf.cluster_idx, nt(1), plot_centroid, nodeSz, mapSz, sigma, wf_cmap);
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(1).Children(1:end-1))
set(nt(1),'XLim',[50 250], 'YLim', [50 250])

% UMAP cluster WFs
nt(2) = nexttile(tl,4,[1 3]);
p = plot(xx_wf,interp_wfs,"LineWidth",1);
nt(2).ColorOrder = wf_cmap;
set(nt(2),'Box','off','XTick',0:1:3,'XLim',[0 2.5],'TickDir','out','YLim',[-1 0.6])
nt(2).YAxis.Visible = 'off';
xlabel("Time [ms]");
leg(1) = legend(nt(2),"N = " + [677,831,167,1147],'Location','southeast','Box','off'); 
leg(1).ItemTokenSize = [4 4]; leg(1).Position([1,2]) = [leg(1).Position(1) + 0.025, leg(1).Position(2) - 0.02];

% WF confusion matrix
nt(3) = nexttile(tl,7,[1 3]);
i = imagesc(norm_cm_wf);
colormap(nt(3),flipud(magma()))
cb(1) = colorbar('Location','eastoutside','Limits',[0 1]); cb(1).Position(3) = 0.01;
set(nt(3),'DataAspectRatio', [2 2 1],'XTick',1:4,'YTick',1:4,'CLim',[0 1],'TickLength',[0 0])
xlabel('Predicted class'); ylabel('True class'); title(' ')
cb(1).Position(1) = sum(nt(3).Position([1,3])) + cb_offset;
cb(1).Title.String = '[%] of True class'; cb(1).Ticks = [0.025 0.95]; cb(1).TickLabels = [0 1];
ticklabels = arrayfun(@(x) sprintf('\\color[rgb]{%f %f %f} %i',wf_cmap(x,1),wf_cmap(x,2),wf_cmap(x,3),x),1:length(nt(3).XTickLabel),'un',0);
nt(3).XTickLabel = ticklabels; nt(3).YTickLabel = ticklabels;

% WF features importance
nt(4) = nexttile(tl,10,[1 3]);
nt(4).ColorOrder = wf_cmap;
p = plot(nt(4),xx_wf,interp_wfs');
arrayfun(@(x) set(x, 'Color', [x.Color 0.3]),p)
hold on
s = scatter(xx_wf,interp_mean_wf,1,yy_wf_pI,'filled');
colormap(nt(4),plasma())
set(nt(4),'Box','off','XTick',0:1:3,'XLim',[0 2.5],'TickDir','out','DataAspectRatio',[1.2 1 1])
nt(4).YAxis.Visible = 'off';
xlabel('Time [ms]')
cb(2) = colorbar('Location','north'); cb(2).Position(4) = 0.01; 
cb(2).Ruler.TickLabelGapOffset = 0; cb(2).Ticks = [0.2 9]; cb(2).TickLabels = [0 9];
text(nt(4),1.3,-0.5,["Average" "waveform"],Color=nt(4).Colormap(1,:),FontSize=fontsz)

% UMAP Combined features
nt(5) = nexttile(tl,13,[1 3]);
RecordingGroup.plot_cluster_outlines(comb.umap, comb.cluster_idx, nt(5), plot_centroid, nodeSz, mapSz, sigma, comb_cmap);
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(5).Children(1:end-1))
set(nt(5),'XLim',[50 250], 'YLim', [50 250])

% Exemplary activity traces
nt(6) = nexttile(tl,16,[1 3]);
scatter(nt(6),scatter_data_comb.spk_t, scatter_data_comb.spk_id, 1, scatter_data_comb.color_array,'filled')
hold on
yline(nt(6),cumsum(scatter_data_comb.clust_sz))
xline(nt(6),230)
ticklabels = arrayfun(@(x) sprintf('\\color[rgb]{%f %f %f} %i',comb_cmap(x,1),comb_cmap(x,2),comb_cmap(x,3),x),1:size(comb_cmap,1),'un',0);
set(nt(6),'XTick',[200 230 260],'XTickLabel',[0 30],'XLim',[200 230],'YTick',scatter_data_comb.ytick_vals,'YTickLabel',ticklabels,'YLim',[0 sum(scatter_data_comb.clust_sz)],'Colormap',comb_cmap)
xlabel('Time [s]'); ylabel('Cluster ID')

% Activity confusion matrix
nt(7) = nexttile(tl,19,[1 3]);
i = imagesc(norm_cm_comb);
colormap(nt(7),flipud(magma()))
cb(4) = colorbar('Location','eastoutside','Limits',[0 1]); cb(4).Position(3) = 0.01;
set(nt(7),'DataAspectRatio', [2 2 1],'XTick',1:size(norm_cm_comb,2),'YTick',1:size(norm_cm_comb,1),'CLim',[0 1],'TickLength',[0 0],'XTickLabelRotation',0)
xlabel('Predicted class'); ylabel('True class'); 
cb(4).Position(1) = sum(nt(7).Position([1,3])) + cb_offset;
cb(4).Ticks = [0.025 0.95]; cb(4).TickLabels = [0 1]; cb(4).Ruler.TickLabelGapOffset = 0;%cb(4).Title.String = '[%] of True Class'; 
ticklabels = arrayfun(@(x) sprintf('\\color[rgb]{%f %f %f} %i',comb_cmap(x,1),comb_cmap(x,2),comb_cmap(x,3),x),1:length(nt(7).XTick),'un',0);
nt(7).XTickLabel = ticklabels; nt(7).YTickLabel = ticklabels; 

% Activity feature heatmap
ax(1) = axes('Position',[0.82    0.62 0.11 0.11]);
imagesc(act.feature_heatmap)
set(ax(1), 'XTick', 1:size(act.feature_heatmap,2), 'YTick', 1:size(act.feature_heatmap,1),'TickLength',[0 0],'DataAspectRatio',[2 2 1],'XTickLabelRotation',0)
yticklabels(act_features); xlabel('Class ID')
colormap(ax(1),plasma())
cb(3) = colorbar(ax(1),'Location','northoutside'); cb(3).Position = [ax(1).Position(1), sum(ax(1).Position([2,4])),  ax(1).Position(3), 0.01];  
cb(3).Ruler.TickLabelGapOffset = -1; cb(3).Title.String = 'Z-score'; cb(3).Title.Position = [24.5 12 0];
set(ax(1),'CLim',[-max(max(abs(act.feature_heatmap))) max(max(abs(act.feature_heatmap)))])

% Activity feature importances
ax(2) = axes('Position',[sum(ax(1).Position([1,3])) + 0.005, ax(1).Position(2)+0.01 0.04 0.09]);
barh(act.pI.sorted(5:-1:1),'FaceColor','k','EdgeColor','none','BarWidth',0.8)
box off; ax(2).YLim = [0.5 size(act.feature_heatmap,1) + 0.5]; ax(2).YTick = [];%ax(2).YAxis.Visible = 'off';
ttl = title('Importance','FontSize',fontsz,'FontWeight','normal','Rotation',-90); ttl.Position = [25 3 0];

%  WF prediction
nt(9) = nexttile(tl,25,[1 3]);
bar(clust_prop_x, norm_wf_pred([1:4,7:10],:),'stacked')
set(nt(9),'Box','off','ColorOrder',wf_cmap,'YTick',[0.05 0.5 1],'YTickLabel',[0 0.5 1])
ylabel('Cluster proportions'); title(nt(9),'Waveform','FontWeight','normal')

% Combined prediction
nt(10) = nexttile(tl,28,[1 3]);
bar(clust_prop_x, norm_comb_pred([1:4,7:10],:),'stacked')
set(nt(10),'Box','off','ColorOrder',comb_cmap)
ylabel(nt(10),''); title(nt(10),'Waveform + activity','FontWeight','normal')
arrayfun(@(x) set(x,'XTickLabel',clust_prop_x_labels,'XTickLabelRotation',0,'TickLength',[0 0],'YLim',[0 1.2]),nt(9:10))
nt(9).XTick = [0 clust_prop_x];
nt(9).XTickLabel = ["LNA    \newlineLine  ",clust_prop_x_labels];
xls = arrayfun(@(x) xline(x, 6.5,'Label','Cell composition','LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom','FontSize',fontsz-1),nt(9:10),'un',0);
arrayfun(@(x) text(x,3,1.18, 'Hetero','VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fontsz-1),nt(9:10));
arrayfun(@(x) text(x,3,1.08, 'geneous','VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fontsz-1),nt(9:10));
arrayfun(@(x) text(x,10,1.18, 'Homo','VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fontsz-1),nt(9:10))
arrayfun(@(x) text(x,10,1.08, 'geneous','VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fontsz-1),nt(9:10))

% QP baseline UMAP + FR change
nt(11) = nexttile(tl,31,[1 3]);
% umap_reduction = RecordingGroup.plot_cluster_outlines(bl_reduction, bl_cluster_idx, nt(11), plot_centroid, nodeSz, mapSz, sigma);
% % arrayfun(@(x) set(nt(11).Children(x),'Visible','off'),1:length(nt(11).Children) - 1)
% % arrayfun(@(x) set(nt(11).Children(x),'MarkerEdgeColor','none'),1:length(nt(11).Children) - 1)
% arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(11).Children(1:end-1))
% % scatter(nt(11), umap_reduction(fr_idx == 1,1),umap_reduction(fr_idx == 1,2),2,repmat([0 0 0],sum(fr_idx == 1),1),'o', 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5)
% scatter(nt(11), umap_reduction(fr_idx == 1,1),umap_reduction(fr_idx == 1,2),1,repmat([1 1 1],sum(fr_idx == 1),1),'.','filled', 'MarkerEdgeColor','k')
% set(nt(11),'XLim',[50 250], 'YLim', [50 250])
% scatter_colors = arrayfun(@(x) x.CData(1,:),nt(11).Children(2:end-1),'un',0);
% scatter_colors = flipud(vertcat(scatter_colors{:}));
% 
% arrayfun(@(x) set(x,'HandleVisibility','off'),nt(11).Children(2:end-2))
% leg(2) = legend(nt(11),["Cluster ID", "FR decrease"],'Box','off','Location','south','NumColumns',2);
% leg(2).ItemTokenSize = [4 4]; leg(2).Position(2) = leg(2).Position(2)-0.06;

% FR change by cluster
nt(12) = nexttile(tl,34,[1 3]);
% ratios = fr_counts./clust_counts * 100;
% bar(ratios)
% box off
% xlabel(nt(12), 'Cluster ID')
% ylabel(nt(12), 'Units with FR decrease [%]')
% title(nt(12), 'FR decrease after QP','FontWeight','normal')
% set(nt(12).Children, 'CData', scatter_colors)
% nt(12).Children.FaceColor = 'flat';
% ticklabels = arrayfun(@(x) sprintf('\\color[rgb]{%f %f %f} %i',scatter_colors(x,1),scatter_colors(x,2),scatter_colors(x,3),x),1:size(scatter_colors,1),'un',0);
% xticklabels(nt(12),ticklabels)
% arrayfun(@(x) text(nt(12), x, 31, sprintf('%i/%i', fr_counts(x), clust_counts(x)),'FontSize',fontsz-1,'Rotation',90),1:length(ratios))

%%
% Finally adjust some positions that are messed up in the process
setAllFontSizes(f5, fontsz)
arrayfun(@(x) set(x,'TickLength',0),cb)
cb(1).Position([2,4]) = nt(3).Position([2,4]);
cb(2).Position(2) = sum(nt(4).Position([2,4])) - 0.01; cb(2).Title.String = 'Predictor importance';
% cb(2).Title.Position(2) =  cb(2).Title.Position(2) - 1;
cb(4).Position([2,4]) = nt(6).Position([2,4]);
cb(4).Position(1) = cb(1).Position(1);
% cb(5).Position([2,4]) = nt(12).Position([2,4]);

%%
exportgraphics(f5,'/home/phornauer/Git/DeePhys/Plots/Figure5_wf_corr.tif','Resolution',300)