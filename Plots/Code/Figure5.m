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
load(fullfile(save_root,'single_cell_predictions.mat'))
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
clust_prop_x = [1:2,4:5,7:8,11:12,14:15];
clust_prop_x_labels = ["     -\newline     WT", "+","      -\newline     SNCA", "+", "     -\newline     GBA", "+", "     -\newline     WT", "+", "      -\newline     A53T", "+"];

%%
f5 = figure('Color','w','Position',[1200 200 595 595],'Units','points');
tl = tiledlayout(5,12,'TileSpacing','tight','Padding','tight');

% UMAP WF clustering
nt(1) = nexttile(tl,1,[1 3]);
RecordingGroup.plot_cluster_outlines(wf.umap, wf.cluster_idx, nt(1), plot_centroid, nodeSz, mapSz, sigma, wf_cmap)
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
leg(1).ItemTokenSize = [4 4]; leg(1).Position([1,2]) = [leg(1).Position(1) + 0.025, leg(1).Position(2) - 0.01];

% WF confusion matrix
nt(3) = nexttile(tl,7,[1 3]);
i = imagesc(norm_cm_wf);
colormap(nt(3),flipud(magma()))
cb(1) = colorbar('Location','eastoutside','Limits',[0 1]); cb(1).Position(3) = 0.01;
set(nt(3),'DataAspectRatio', [2 2 1],'XTick',1:4,'YTick',1:4,'CLim',[0 1],'TickLength',[0 0])
xlabel('Predicted Class'); ylabel('True Class'); title(' ')
cb(1).Position(1) = sum(nt(3).Position([1,3])) + cb_offset;
cb(1).Title.String = '[%] of True Class'; cb(1).Ticks = [0.025 0.95]; cb(1).TickLabels = [0 1];
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
text(nt(4),1.3,-0.25,["Average" "Waveform"],Color=nt(4).Colormap(1,:),FontSize=fontsz)

% UMAP activity features
nt(5) = nexttile(tl,13,[1 3]);
RecordingGroup.plot_cluster_outlines(act.umap, act.cluster_idx, nt(5), plot_centroid, nodeSz, mapSz, sigma, act_cmap)
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(5).Children(1:end-1))
set(nt(5),'XLim',[50 250], 'YLim', [50 250])

% Exemplary activity traces
nt(6) = nexttile(tl,16,[1 3]);
scatter(scatter_data.spk_t, scatter_data.spk_id, 1, scatter_data.color_array,'filled')
hold on
yline(nt(6),cumsum(scatter_data.clust_sz))
xline(nt(6),230)
ticklabels = arrayfun(@(x) sprintf('\\color[rgb]{%f %f %f} %i',act_cmap(x,1),act_cmap(x,2),act_cmap(x,3),x),1:size(act_cmap,1),'un',0);
set(nt(6),'XTick',[200 230 260],'XTickLabel',[0 30 60],'XLim',[200 230],'YTick',scatter_data.ytick_vals,'YTickLabel',ticklabels,'YLim',[0 sum(scatter_data.clust_sz)],'Colormap',act_cmap)
xlabel('Time [s]'); ylabel('Cluster ID')

% Activity confusion matrix
nt(7) = nexttile(tl,19,[1 3]);
i = imagesc(norm_cm_act);
colormap(nt(7),flipud(magma()))
cb(4) = colorbar('Location','eastoutside','Limits',[0 1]); cb(4).Position(3) = 0.01;
set(nt(7),'DataAspectRatio', [2 2 1],'XTick',1:size(norm_cm_act,2),'YTick',1:size(norm_cm_act,1),'CLim',[0 1],'TickLength',[0 0],'XTickLabelRotation',0)
xlabel('Predicted Class'); ylabel('True Class'); 
cb(4).Position(1) = sum(nt(7).Position([1,3])) + cb_offset;
cb(4).Ticks = [0.025 0.95]; cb(4).TickLabels = [0 1]; %cb(4).Title.String = '[%] of True Class'; 
ticklabels = arrayfun(@(x) sprintf('\\color[rgb]{%f %f %f} %i',act_cmap(x,1),act_cmap(x,2),act_cmap(x,3),x),1:length(nt(7).XTick),'un',0);
nt(7).XTickLabel = ticklabels; nt(7).YTickLabel = ticklabels;

% Activity feature heatmap
ax(1) = axes('Position',[0.82    0.62 0.11 0.11]);
imagesc(act.feature_heatmap)
set(ax(1), 'XTick', 1:size(act.feature_heatmap,2), 'YTick', 1:size(act.feature_heatmap,1),'TickLength',[0 0],'DataAspectRatio',[2 2 1],'XTickLabelRotation',0)
yticklabels(act_features); xlabel('Class ID')
colormap(ax(1),plasma())
cb(3) = colorbar(ax(1),'Location','northoutside'); cb(3).Position = [ax(1).Position(1), sum(ax(1).Position([2,4])),  ax(1).Position(3), 0.01];  
cb(3).Ruler.TickLabelGapOffset = -1; cb(3).Title.String = 'Normalized value'; cb(3).Title.Position = [24.5 12 0];
set(ax(1),'CLim',[-max(max(abs(act.feature_heatmap))) max(max(abs(act.feature_heatmap)))])

% Activity feature importances
ax(2) = axes('Position',[sum(ax(1).Position([1,3])) + 0.005, ax(1).Position(2)+0.01 0.04 0.09]);
barh(act.pI.sorted(5:-1:1),'FaceColor','k','EdgeColor','none','BarWidth',0.8)
box off; ax(2).YLim = [0.5 size(act.feature_heatmap,1) + 0.5]; ax(2).YTick = [];%ax(2).YAxis.Visible = 'off';
ttl = title('Importance','FontSize',fontsz,'FontWeight','normal','Rotation',-90); ttl.Position = [25 3 0];

% Combined UMAP plot
nt(9) = nexttile(tl,25,[1 3]);
RecordingGroup.plot_cluster_outlines(comb.umap, comb.cluster_idx, nt(9), plot_centroid, nodeSz, mapSz, sigma, comb_cmap)
arrayfun(@(x) set(x,'MarkerEdgeAlpha',0.2),nt(9).Children(1:end-1))
set(nt(9),'XLim',[50 250], 'YLim', [50 250])

% Activity combined clusters
nt(10) = nexttile(tl,28,[1 3]);
scatter(nt(10),scatter_data_comb.spk_t, scatter_data_comb.spk_id, 1, scatter_data_comb.color_array,'filled')
hold on
yline(nt(10),cumsum(scatter_data_comb.clust_sz))
xline(nt(10),230)
ticklabels = arrayfun(@(x) sprintf('\\color[rgb]{%f %f %f} %i',comb_cmap(x,1),comb_cmap(x,2),comb_cmap(x,3),x),1:size(comb_cmap,1),'un',0);
set(nt(10),'XTick',[200 230 260],'XTickLabel',[0 30],'XLim',[200 230],'YTick',scatter_data_comb.ytick_vals,'YTickLabel',ticklabels,'YLim',[0 sum(scatter_data_comb.clust_sz)],'Colormap',comb_cmap)
xlabel('Time [s]'); ylabel('Cluster ID')

% Waveforms
nt(11) = nexttile(tl,34,[1 3]);
p = plot(nt(11),xx_wf,interp_comb_wfs,"LineWidth",1);
nt(11).ColorOrder = comb_cmap;
set(nt(11),'Box','off','XTick',0:1:3,'XLim',[0 2.5],'TickDir','out','YLim',[-1 0.6],'DataAspectRatio',[1 1 1])
nt(11).YAxis.Visible = 'off';
xlabel("Time [ms]");
% leg(1) = legend(nt(2),"N = " + [1230,1118,1107,392],'Location','southeast','Box','off'); 
% leg(1).ItemTokenSize = [4 4]; leg(1).Position([1,2]) = [leg(1).Position(1) + 0.02, leg(1).Position(2) - 0.02];

% Confusion matrix
nt(12) = nexttile(tl,31,[1 3]);
i = imagesc(norm_cm_comb);
colormap(nt(12),flipud(magma()))
cb(5) = colorbar('Location','eastoutside','Limits',[0 1]); cb(5).Position(3) = 0.01;
set(nt(12),'DataAspectRatio', [2 2 1],'XTick',1:size(norm_cm_comb,2),'YTick',1:size(norm_cm_comb,1),'CLim',[0 1],'TickLength',[0 0],'XTickLabelRotation',0)
xlabel('Predicted Class'); ylabel('True Class'); 
cb(5).Position(1) = sum(nt(12).Position([1,3])) + cb_offset;
cb(5).Ticks = [0.025 0.95]; cb(5).TickLabels = [0 1]; cb(5).Ruler.TickLabelGapOffset = 0;%cb(4).Title.String = '[%] of True Class'; 
ticklabels = arrayfun(@(x) sprintf('\\color[rgb]{%f %f %f} %i',comb_cmap(x,1),comb_cmap(x,2),comb_cmap(x,3),x),1:length(nt(12).XTick),'un',0);
nt(12).XTickLabel = ticklabels; nt(12).YTickLabel = ticklabels; 

% WF prediction
nt(13) = nexttile(tl,37,[1 4]);
bar(clust_prop_x, norm_wf_pred,'stacked')
set(nt(13),'Box','off','ColorOrder',wf_cmap,'YTick',[0.05 0.5 1],'YTickLabel',[0 0.5 1])
ylabel('Cluster Proportions'); title(nt(13),'Waveform Classification','FontWeight','normal')

% Activity prediction
nt(14) = nexttile(tl,41,[1 4]);
bar(clust_prop_x, norm_act_pred,'stacked')
set(nt(14),'Box','off','ColorOrder',act_cmap)
ylabel(''); title(nt(14),'Activity Classification','FontWeight','normal')

% Combined prediction
nt(15) = nexttile(tl,45,[1 4]);
bar(clust_prop_x, norm_comb_pred,'stacked')
set(nt(15),'Box','off','ColorOrder',comb_cmap)
ylabel(nt(15),''); title(nt(15),'Combined Classification','FontWeight','normal')
arrayfun(@(x) set(x,'XTickLabel',clust_prop_x_labels,'XTickLabelRotation',0,'TickLength',[0 0],'YLim',[0 1.2]),nt(13:15))
nt(13).XTick = [0 clust_prop_x];
nt(13).XTickLabel = ["LNA    \newlineLine  ",clust_prop_x_labels];
xls = arrayfun(@(x) xline(x, 9.5,'Label','Cell Composition','LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom','FontSize',fontsz),nt(13:15),'un',0);
arrayfun(@(x) text(x,4.5,1.1, 'Heterogeneous','VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fontsz),nt(13:15))
arrayfun(@(x) text(x,13,1.1, 'Homogeneous','VerticalAlignment','middle','HorizontalAlignment','center','FontSize',fontsz),nt(13:15))

%%
% Finally adjust some positions that are messed up in the process
setAllFontSizes(f5, fontsz)
arrayfun(@(x) set(x,'TickLength',0),cb)
cb(1).Position([2,4]) = nt(3).Position([2,4]);
cb(2).Position(2) = sum(nt(4).Position([2,4])) - 0.01; cb(2).Title.String = 'Feature Importance';
% cb(2).Title.Position(2) =  cb(2).Title.Position(2) - 1;
cb(4).Position([2,4]) = nt(6).Position([2,4]);
cb(5).Position([2,4]) = nt(12).Position([2,4]);

ax(3) = axes('Position',[0.65 0.43 0.07 0.06]);
p = plot(ax(3),xx_wf,interp_comb_wfs,"LineWidth",1);
set(ax(3),'XTick',[],'YTick',[],'Box','on','XLim',[0.75 0.9],'YLim',[-0.1 0.55],'ColorOrder',comb_cmap,'Position',[0.55 + 0.25 0.435 0.06 0.06],'XColor',[0.5 0.5 0.5], 'YColor',[0.5 0.5 0.5])
r = rectangle(nt(11),'Position',[0.75 -0.1 0.15 0.65],'EdgeColor',[0.5 0.5 0.5]);

ax(4) = axes;
p = plot(ax(4),xx_wf,interp_comb_wfs,"LineWidth",1);
set(ax(4),'XTick',[],'YTick',[],'Box','on','XLim',[0.95 1.2],'YLim',[-0.6 0.2],'ColorOrder',comb_cmap,'Position',[0.65 + 0.25 0.435 0.08 0.06],'XColor',[0.5 0.5 0.5], 'YColor',[0.5 0.5 0.5])
r = rectangle(nt(11),'Position',[0.95 -0.6 0.25 0.85],'EdgeColor',[0.5 0.5 0.5]);

%%
exportgraphics(f5,'/home/phornauer/Git/DeePhys/Plots/Figure5.tif','Resolution',300)

%%
figure('Color','w','Position',[800 800 200 250]);
imagesc(full_test_acc)
xticks(1:size(full_test_acc,2))
xticklabels(["3","4","5","All"])
yticks(1:size(full_test_acc,1))
yticklabels(feature_groups)
xlabel('Week')
set(gca,'TickLength',[0 0])
colormap(greens)
colorbar('Limits',[0 1])
title('Classification accuracy','FontWeight','normal')
setAllFontSizes(gcf,6.5)