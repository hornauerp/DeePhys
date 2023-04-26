addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Set some parameters
greens = othercolor('Greens9',100); %greens = greens(30:end,:);
heatmap_colormap = othercolor('RdBu9',40); heatmap_colormap((size(heatmap_colormap,1)/2)-2:size(heatmap_colormap,1)/2+2,:) = [];
colors_asyn = [0 0 0; 0.5 0.5 0.5; 1 1 1];
fontsz = 6.5;
capsz = 3;
boxsize = 0.0177; 
x_loc = 0.1; %Offset for first plot 
y_loc = 0.2949;  %Offset 
x_plot = 0.12; %Size for timeline plots
y_plot = 0.08; %Size for timeline plots
diff_heatmaps = 0.01; %
diff_colorbar = 0.01;
cb_height = 0.01;
TAKE_N_BEST = 10;
tps = 6:7:27; %Timepoints
coactivity_lim = [0 700];
% labels = [repelem("TEST",10), "", repelem("TEST",10)];
labels = ["sRF","sRM","sEAF","sAMI","sSFR","CVI","sLPF","sMEF","sPAM","sEFD","",... %SC
    "nPMW","nCCD","nRF","nSFR","nSFD","nMD5","nMD10","nFMA","nTCT","nFMI"]; %NW
example_idx = [11 + 3, 11 + 1,6,3]; %Descending for right plotting order
example_labels = ["nRF [Hz]","nPMW","CVI","sEAF"];

%% Load data
save_root = "/home/phornauer/Git/DeePhys/Data/Figure4";
load(fullfile(save_root,'lna_accuracy.mat'))
load(fullfile(save_root,'FRET_results.mat'))
load(fullfile(save_root,'lna_comparisons.mat'))
load(fullfile(save_root,'lna_sc_umap.mat'))
load(fullfile(save_root,'lna_nw_umap.mat'))
load(fullfile(save_root,'lna_comb_umap.mat'))
load(fullfile(save_root,'coactivity_plot.mat'))
load(fullfile(save_root,'lna_age_prediction.mat'))

%% Prepare some data
p_vals_wt = min([wt_p_G;wt_p_GxT]);
wt_discretized = 4 - discretize(p_vals_wt,[0,0.001, 0.01, 0.05 Inf]);
wt_stars = [wt_discretized(1:10), 0, wt_discretized(11:20)];

p_vals_a53t = min([a53t_p_G;a53t_p_GxT]);
a53t_discretized = 4 - discretize(p_vals_a53t,[0,0.001, 0.01, 0.05 Inf]);
a53t_stars = [a53t_discretized(1:10), 0, a53t_discretized(11:20)];

timeline_data = permute(cat(1,means(1:10,:,:),zeros(1,4,4),means(11:20,:,:)), [3 2 1]);
timeline_sd = permute(cat(1,sds(1:10,:,:),zeros(1,4,4),sds(11:20,:,:)), [3 2 1]);
jitter = linspace(-1,1,size(timeline_data,1)); 

wt_log_ratio = [wt_color_mat(1:10,:); zeros(1,size(wt_color_mat, 2)); wt_color_mat(11:20,:)]';
a53t_log_ratio = [a53t_color_mat(1:10,:); zeros(1,size(a53t_color_mat, 2)); a53t_color_mat(11:20,:)]';

%% Main figure
f4 = figure('Color','w','Position',[1200 200 595 595],'Units','points');
% tl = tiledlayout(4,10,'TileSpacing','compact','Padding','compact');
% Leave first tiles free for staining images

% asyn FRET bar plot
% nt(1) = nexttile(9); nt(1).Layout.TileSpan = [1 2];
syn_bar = axes('Position',[0.8249    0.8    0.1276    0.17]);
bars = bar([norm_means(1:3);norm_means(4:6)],'BarWidth',0.7);
arrayfun(@(x) set(bars(x),'FaceColor', colors_asyn(x,:)),1:length(bars));
hold on
e = errorbar([bars.XEndPoints],[bars.YData],norm_sds(:),'k','linestyle','none','CapSize',capsz);
s = plot(syn_bar, repmat(bars(end).XEndPoints,[3 1]), e.YData(5:6) + e.YPositiveDelta(5:6) + repmat(0.06:0.06:0.18,[2 1])','k*','MarkerSize',3);

box off
xticklabels(["WT","A53T"])
ylabel(['Normalized ' char(945) '-syn [%]'])

ylim([0 1.50]); yticks(sort([get(gca,'YLim'), 1])); yticklabels(get(gca,'YTick') * 100)
leg_asyn = legend(["Untreated","ntLNA","LNA"],'Box','off','Location','north'); leg_asyn.ItemTokenSize = [4 4]; leg_asyn.Position(2) = 0.925;
% title(['Total ', char(945), '-syn'],'FontWeight','normal','FontSize',fontsz)

sp = stackedplot(spike_tbl,'XLabel',"Time [s]",'XLimits',coactivity_lim, 'Position',[0.0830    0.29    0.1204    0.415]);
arrayfun(@(x) set(sp.LineProperties(x),'Color',stacked_color(x,:)),1:length(sp.LineProperties));
ax = findobj(sp.NodeChildren, 'Type','Axes');

% Break here because it gives an error otherwise
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','FontSize',fontsz)
arrayfun(@(x) set(ax(x).YLabel, 'String', y_labels(x)),1:length(ax))
% 
arrayfun(@(x) set(x, 'XTick',0:300:600, 'FontSize', fontsz, 'XTickLabels',0:30:60,'YLim',[-100 500],'YTick',[0 400],'YTickLabels',[0 400],'YTickLabelRotation',0),ax)
arrayfun(@(x) set(sp.AxesProperties(x),'YLimits',[-50 500]),1:length(sp.LineProperties));
ax(1).XLim = [0 700];
set(ax(4).Title,'String',"Network coactivity" + newline + "[spikes/100ms]",'FontWeight','normal','FontSize',fontsz)


wt_hmap = axes('Position',[0.3,...
                    0.2949,...
                    size(wt_log_ratio,1)*boxsize,...
                    size(wt_log_ratio,2)*boxsize]);
hmap = imagesc(wt_log_ratio');
hold on
%Mark plotted features in heatmap
arrayfun(@(x) rectangle('Position',[0.5 x-0.5 size(wt_log_ratio,1) 1],'FaceColor','none','LineWidth',1.5),example_idx);
wt_hmap.TickLength = [0 0];
wt_hmap.FontSize = fontsz;
wt_hmap.YTick = 1:length(labels);
wt_hmap.YTickLabel = labels;
wt_hmap.XTick = 1:5;
wt_hmap.XLabel.String = "Week";
xtickangle(wt_hmap,0)
title("WT")

colormap(heatmap_colormap)
set(gca,'CLim',[log(1/3) log(3)])
cb = colorbar('Position',[wt_hmap.Position(1) + wt_hmap.Position(3)*0.5, wt_hmap.Position(2)+wt_hmap.Position(4)+0.03, wt_hmap.Position(3)*2, cb_height],...
    'Location','northoutside');
cb.Ticks = [-1.05 0 1.05];
cb.TickLabels = ["<0.3", "1", ">3"];
cb.Ruler.TickLabelGapOffset = -1;
cb.Ruler.TickLabelRotation = 0;
cb.Title.String = "\bf Control/LNA";
cb.Title.Position = [30 14 0];

wt_ss = axes('Position',[wt_hmap.Position(1)+wt_hmap.Position(3)+0.015,...
    y_loc,...
    0.05,...
    wt_hmap.Position(4)]);
plot(ones(1,sum(wt_stars>0)),length(labels)-find(wt_stars)+1,'k*','MarkerSize',3)
hold on
plot(ones(1,sum(wt_stars>1))+0.1,length(labels)-find(wt_stars>1)+1,'k*','MarkerSize',3)
plot(ones(1,sum(wt_stars>2))+0.2,length(labels)-find(wt_stars>2)+1,'k*','MarkerSize',3)
ylim([0.5 length(labels)+0.5])
xlim([1 1.5])
t_MLM = title('LMM','FontWeight','normal','FontSize',fontsz,'HorizontalAlignment','right');
t_MLM.Position(1) = 1.27;
axis off
annotation(f4,'rectangle',...
    [wt_hmap.Position(1)+wt_hmap.Position(3)+0.005    wt_hmap.Position(2)    0.0400    wt_hmap.Position(4)],...
    'LineStyle',':','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.1);

% Draw gap between feature groups
feature_sep1 = yline(wt_hmap,TAKE_N_BEST+0.5,'LineWidth',1);
feature_sep2 = yline(wt_hmap,TAKE_N_BEST+1.5,'LineWidth',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% A53T %%%%%%%%%%%%%%%%%%%%%%%%%%%%
a53t_hmap = axes('Position',[wt_ss.Position(1)+ 0.05,...
                    0.2949,...
                    size(a53t_log_ratio,1)*boxsize,...
                    size(a53t_log_ratio,2)*boxsize]);
hmap = imagesc(a53t_log_ratio');
hold on
%Mark plotted features in heatmap
arrayfun(@(x) rectangle('Position',[0.5 x-0.5 size(wt_log_ratio,1) 1],'FaceColor','none','LineWidth',1.5),example_idx);
a53t_hmap.TickLength = [0 0];
a53t_hmap.FontSize = fontsz;
a53t_hmap.YTick = 1:length(labels);
a53t_hmap.YTickLabel = [];
a53t_hmap.XTick = 1:5;
a53t_hmap.XLabel.String = "Week";
xtickangle(a53t_hmap,0)
title("A53T")

colormap(heatmap_colormap)
set(gca,'CLim',[log(1/3) log(3)])

a53t_ss = axes('Position',[a53t_hmap.Position(1)+a53t_hmap.Position(3)+0.015,...
    y_loc,...
    0.05,...
    a53t_hmap.Position(4)]);
plot(ones(1,sum(a53t_stars>0)),length(labels)-find(a53t_stars)+1,'k*','MarkerSize',3)
hold on
plot(ones(1,sum(a53t_stars>1))+0.1,length(labels)-find(a53t_stars>1)+1,'k*','MarkerSize',3)
plot(ones(1,sum(a53t_stars>2))+0.2,length(labels)-find(a53t_stars>2)+1,'k*','MarkerSize',3)
ylim([0.5 length(labels)+0.5])
xlim([1 1.5])
t_MLM = title('LMM','FontWeight','normal','FontSize',fontsz,'HorizontalAlignment','right');
t_MLM.Position(1) = 1.27;
axis off
annotation(f4,'rectangle',...
    [a53t_hmap.Position(1) + a53t_hmap.Position(3) + 0.005,   a53t_hmap.Position(2),    0.0400,    a53t_hmap.Position(4)],...
    'LineStyle',':','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.1);

% Draw gap between feature groups
white_gap = annotation(f4,'rectangle',...
    [wt_hmap.Position(1) 0.473 0.255 0.016],...
    'LineStyle','none',...
    'FaceColor',[1 1 1]);
feature_sep1 = yline(a53t_hmap,TAKE_N_BEST+0.5,'LineWidth',1);
feature_sep2 = yline(a53t_hmap,TAKE_N_BEST+1.5,'LineWidth',1);


for p = 1:length(example_idx)
    exp_plots{p} = axes('Position',...
        [white_gap.Position(1)+white_gap.Position(3)+0.08,...
        y_loc+(p-1)*(y_plot+0.03),...
        x_plot,...
        y_plot]);
    for i = 1:size(timeline_data,1)
        x = tps+1+jitter(i);
        y = timeline_data(i,:,example_idx(p));
        xx = linspace(min(x),max(x),100); yy = pchip(x,y,xx);
        plot(xx,yy,'Color',stacked_color(i,:))
        hold on
        errorbar(tps+1+jitter(i),timeline_data(i,:,example_idx(p)),timeline_sd(i,:,example_idx(p)),...
            'LineWidth',1,'Color',stacked_color(i,:),'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',2,...
            'MarkerFaceColor',stacked_color(i,:),'HandleVisibility','off');
        set(gca,'FontSize',fontsz)
        marg = get(gca,'ylabel');
        set(marg,'Margin',3)
    end
    box off
    if p == 1
        xlabel('Week')
        xticklabels(1:4)
    else
        exp_plots{p}.XTickLabel = [];
    end
    if p == 4
        leg = legend(["WT", "WT + LNA", "A53T", "A53T + LNA"],'box','off');
        leg.ItemTokenSize = [5 5];
        leg.NumColumns = 2;
        leg.Position(1) = exp_plots{p}.Position(1)-0.055;
        leg.Position(2) = exp_plots{p}.Position(2) + exp_plots{p}.Position(4) + 0.005;
%         leg.Position = [0.52, 0.9, 0.065, 0.0445];
    end
    ylabel(example_labels(p))
    xlim([3 max(xlim)])
    xticks(7:7:28)
    
    yt = get(gca,'YTick');
    yticks(linspace(min(yt),max(yt),2)) %Change here to insert more ticks
    yl = get(gca,'YLim');
    ylim([yl(1)-max(abs(yl))*0.1 yl(2)])
end

plot_centroid = false;
nodeSz = 10;
mapSz = 300;
sigma = mapSz/60;
cmap = stacked_color([4,3,2,1],:);
y_loc_umap = 0.26;
x_sz = 0.13;
y_sz = 0.13;
umap_sc = axes('Position',...
    [exp_plots{1}.Position(1)+exp_plots{1}.Position(3)+0.05,...
    y_loc_umap+2*(y_sz+0.04),...
    x_sz,...
    y_sz]);

RecordingGroup.plot_cluster_outlines(sc_reduction, sc_true_idx, umap_sc, plot_centroid, nodeSz, mapSz, sigma, cmap)
title('Single-cell features','FontWeight','normal')

ump_nw = axes('Position',...
    [exp_plots{1}.Position(1)+exp_plots{1}.Position(3)+0.05,...
    y_loc_umap+(y_sz+0.04),...
    x_sz,...
    y_sz]);

RecordingGroup.plot_cluster_outlines(nw_reduction, nw_true_idx, ump_nw, plot_centroid, nodeSz, mapSz, sigma, cmap)
title('Network features','FontWeight','normal')

umap_comb = axes('Position',...
    [exp_plots{1}.Position(1)+exp_plots{1}.Position(3)+0.05,...
    y_loc_umap,...
    x_sz,...
    y_sz]);

RecordingGroup.plot_cluster_outlines(comb_reduction, comb_true_idx, umap_comb, plot_centroid, nodeSz, mapSz, sigma, cmap)
title('Combined features','FontWeight','normal')

setAllFontSizes(gcf,fontsz)
exportgraphics(f4,'/home/phornauer/Git/DeePhys/Plots/Figure4a.tif','Resolution',300)

%% Second part
acc_boxsz = 0.03;
y_pred_wt = vertcat(wt_age.Y_pred);
y_test_wt = vertcat(wt_age.Y_test);
y_pred_a53t = vertcat(a53t_age.Y_pred);
y_test_a53t = vertcat(a53t_age.Y_test);
y_pred = [y_pred_wt; y_pred_a53t];
y_test = [y_test_wt; y_test_a53t];
age_color_idx = [zeros(size(y_pred_wt)); ones(size(y_pred_a53t))];
min_age_diff = min(diff(unique(y_test)));
age_x_vals = y_test/min_age_diff;

%%
f4b = figure('Color','w','Position',[1200 200 595 595],'Units','points');
acc_heatmap = axes('Position',[0.1 0.07 size(lna_accuracy_mat,1)*acc_boxsz,...
                    size(lna_accuracy_mat,2)*acc_boxsz]); 
imagesc(acc_heatmap,lna_accuracy_mat(:,[5:8,1:4,end]).Variables')
hold on
yline([8.5 4.5],'LineWidth',1,'Color','k')
xticks(1:size(lna_accuracy_mat,1))
row1 = {'WT','A53T','WT','A53T'}; row2 = {'','','LNA','LNA'}; labelArray = [row1; row2];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
xticklabels(tickLabels)
xticklabels(["WT","A53T","WT LNA","A53T LNA"])
set(acc_heatmap,'XTickLabelRotation',45)
yticks(1:size(lna_accuracy_mat,2))
yticklabels(feature_groups([5:8,1:4,end]))
set(gca,'TickLength',[0 0])
colormap(acc_heatmap, greens)
cb = colorbar('Limits',[0 1],'Position',[0.2300    acc_heatmap.Position(2)    0.0100    0.2689]);
set(cb.Title,"String", "Classification accuracy","Position",[20    60     0],"Rotation", -90);
% title('Classification accuracy','FontWeight','normal')


full_umap = axes('Position',[acc_heatmap.Position(1) + acc_heatmap.Position(3) + 0.11, acc_heatmap.Position(2), 0.25, acc_heatmap.Position(4)]);
RecordingGroup.plot_cluster_outlines(full_lna_reduction, full_lna_idx, full_umap, plot_centroid, nodeSz, mapSz, sigma, cmap)
% title('Best features','FontWeight','normal','FontSize',fontsz)
leg = legend(["A53T + LNA", "A53T", "WT + LNA","WT"],'box','off');
leg.Location = 'northoutside';
leg.ItemTokenSize = [4 4];
leg.NumColumns = 2;
% leg.Position(1) = exp_plots{p}.Position(1)-0.055;
% leg.Position(2) = exp_plots{p}.Position(2) + exp_plots{p}.Position(4) + 0.005;


age_pred = axes('Position',[full_umap.Position(1) + full_umap.Position(3) + 0.07, full_umap.Position(2), 0.25, full_umap.Position(4)]);
bc = boxchart(age_pred, age_x_vals, y_pred,'GroupByColor',age_color_idx,'MarkerSize',1);
hold on
p = plot(age_pred, [0 max(age_x_vals)*1.2],[0 max(y_test)*1.2],'k--'); p.Color(4) = 0.3;
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xticks(1:length(tps))
xlabel('Culture age [weeks]')
yticks(tps)
yticklabels(1:length(tps))
ylabel('Predicted age [weeks]')
set(age_pred,'FontSize',fontsz)
leg = legend("WT + LNA", "A53T + LNA",'Box','off','Location','northwest'); leg.ItemTokenSize = [5 5];
setAllFontSizes(gcf,6.5)

exportgraphics(f4b,'/home/phornauer/Git/DeePhys/Plots/Figure4b.tif','Resolution',300)

