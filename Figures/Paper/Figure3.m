rg_params.Selection.Inclusion = {{'Mutation',"WT"}};
rg_WT = RecordingGroup(fcdi_array,rg_params);
rg_params.Selection.Inclusion = {{'Mutation',"A53T"}};
rg_A53T = RecordingGroup(fcdi_array,rg_params);

%% Scatter plot
time_cutout = [0 120];
f3a = figure('Color','w','Position',[1200 100 595 595],'Units','points');
tl = tiledlayout(7,3,'TileSpacing','compact','Padding','tight');
nt(1) = nexttile(tl,3);
rg_WT.Recordings(end).PlotNetworkScatter(time_cutout,blue)
nt(1).XLabel = [];%title('WT','Color',blue)

nt(2) = nexttile(tl,6);
rg_A53T.Recordings(end).PlotNetworkScatter(time_cutout,red)
%title('A53T','Color',red)
setAllFontSizes(f3a, fontsz)
exportgraphics(f3a,'/home/phornauer/Git/DeePhys/Plots/Figure3a.tif','Resolution',300)

%% SET PARAMETERS
TAKE_N_BEST = 15; %Number of features with highest accuracy to include
N_PRED = 5; %Number of feature importances for age prediction to plot
boxsize = 0.02; 
x_loc = 0.1; %Offset for first plot 
y_loc = 0.3; %Offset 
x_plot = 0.15; %Size for timeline plots
y_plot = 0.1; %Size for timeline plots
diff_heatmaps = 0.01; %
diff_colorbar = 0.01;
cb_height = 0.01;
fontsz = 6.5;
capsz = 3;
tps = 6:7:34; %Timepoints

example_idx = [TAKE_N_BEST + 1 + 4,TAKE_N_BEST + 1 + 3,6,2]; %Descending for right plotting order
example_labels = ["nSFR","nRF [Hz]","CVI","sRM"];
N_groups = 2;

%UMAP plot 
x_sz = 0.2; 
y_sz = 0.2*3/4;

reds = othercolor('Reds9',9);
red = reds(7,:);
greens = othercolor('Greens9',100); greens = greens(30:end,:);
green = greens(7,:);
blues = othercolor('Blues9',9);
blue = blues(7,:);
greys = othercolor('Greys9',20);
greys = greys(6:15,:);

heatmap_colormap = othercolor('RdBu9',N_groups*20); heatmap_colormap((size(heatmap_colormap,1)/2)-2:size(heatmap_colormap,1)/2+2,:) = [];

%% LOAD DATA
save_path = "/home/phornauer/Git/DeePhys/Data/Figure3";
load(fullfile(save_path,'network_p_vals.mat'))
load(fullfile(save_path,'network_heatmap.mat'))
load(fullfile(save_path,'network_values.mat'))
load(fullfile(save_path,'network_features.mat'))
nw_features = feature_names{1};

load(fullfile(save_path,'single_cell_p_vals.mat'))
load(fullfile(save_path,'single_cell_heatmap.mat'))
load(fullfile(save_path,'single_cell_values.mat'))
load(fullfile(save_path,'single_cell_features.mat'))
sc_features = feature_names{1};

load(fullfile(save_path,'single_cell_clustering.mat'))
load(fullfile(save_path,'network_clustering.mat'))
load(fullfile(save_path,'comb_clustering.mat'))

load(fullfile(save_path, 'feature_group_heatmap.mat'))%% Age regression
level = "Recording"; %Unit or Recording
alg = "rf"; %rf, svm, cnb, knn
stratification_var = "Mutation"; %Specify the variable by which to split training and test dataset 
stratification_values = []; %Corresponding to stratification_var, if a value is specified then this will be used as training data (e.g. train on untreated, test on treated)
pooling_vals = {};
network_features = ["all"];
unit_features = ["all"];
useClustered = false;
normalization_var = "PlatingDate";
N_hyper = 0; %If >0 do hyperparameter optimization
K_fold = -1; % number of K-fold CV

%% Age regression WT
rg_params.Selection.Inclusion = {{'Mutation','WT'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',12},{'PlatingDate',200121}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt_rg = RecordingGroup(rec_array, rg_params);

wt_age_result = wt_rg.predictAge(level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold);

load(fullfile(save_path, 'base_age_prediction.mat'))

%% DATA SELECTION
[nw_acc_sorted, nw_sorted_idx] = sort(nw_test_acc,'descend');
nw_sorted_idx = nw_sorted_idx(1:TAKE_N_BEST);
nw_acc_sorted = nw_acc_sorted(1:TAKE_N_BEST);
nw_features_sorted = nw_features(nw_sorted_idx);
nw_pred_imp_sorted = nw_pred_imp(nw_sorted_idx,:);
nw_color_mat_sorted = nw_color_mat(nw_sorted_idx,:);
nw_p_val = min([nw_p_G; nw_p_GxT],[],1);
nw_p_val_sorted = 4 - discretize(nw_p_val(nw_sorted_idx),[0,0.001, 0.01, 0.05 Inf]);
nw_mean_sorted = nw_mean_mat(nw_sorted_idx,:,:);
nw_sd_sorted = nw_sd_mat(nw_sorted_idx,:,:);

[sc_acc_sorted, sc_sorted_idx] = sort(sc_test_acc,'descend');
sc_sorted_idx = sc_sorted_idx(1:TAKE_N_BEST);
sc_acc_sorted = sc_acc_sorted(1:TAKE_N_BEST);
sc_features_sorted = sc_features(sc_sorted_idx);
sc_pred_imp_sorted = sc_pred_imp(sc_sorted_idx,:);
sc_color_mat_sorted = sc_color_mat(sc_sorted_idx,:);
sc_p_val = min([sc_p_G; sc_p_GxT],[],1);
sc_p_val_sorted = 4 - discretize(sc_p_val(sc_sorted_idx),[0,0.001, 0.01, 0.05 Inf]);
sc_mean_sorted = sc_mean_mat(sc_sorted_idx,:,:);
sc_sd_sorted = sc_sd_mat(sc_sorted_idx,:,:);

%% PREPARE INPUT DATA
heatmap_log_ratio = [sc_color_mat_sorted; zeros(1,size(sc_color_mat, 2)); nw_color_mat_sorted]';
labels = string([sc_features_sorted "" nw_features_sorted]);
labels = ["sRF","sRM","sEAF","sAMI","sSFR","CVI","sLPF","sMEF","sPAM","sEFD","PAF","sCCD","sTCT","sSES","sCFS","",... %SC
    "nPMW","nCCD","nRF","nSFR","nSFD","nMD5","nMD10","nFMA","nTCT","nFMI","MIB","nRM","nAMI","nCFS","MFT"]; %NW
labels(example_idx) = "\bf " +labels(example_idx); %Make example features bold
feature_importances = [sc_pred_imp_sorted;zeros(1,size(sc_color_mat,2)); nw_pred_imp_sorted];
feature_importances(feature_importances<=0) = nan; %Convert for scatter plot
accuracy_values = [sc_acc_sorted; 0; nw_acc_sorted'];
num_stars = [sc_p_val_sorted, 0, nw_p_val_sorted]; %Find number of significance stars to plot
jitter = linspace(-1,1,size(heatmap_log_ratio,1)); %Jitter for timeline plots
timeline_data = permute(cat(1,sc_mean_sorted,zeros(1,5,2),nw_mean_sorted), [3 2 1]);
timeline_sd = permute(cat(1,sc_sd_sorted,zeros(1,5,2),nw_sd_sorted),[3 2 1]);

y_pred_wt = vertcat(wt_age.Y_pred);
y_test_wt = vertcat(wt_age.Y_test);
y_pred_a53t = vertcat(a53t_age.Y_pred);
y_test_a53t = vertcat(a53t_age.Y_test);
y_pred = [y_pred_wt; y_pred_a53t];
y_test = [y_test_wt; y_test_a53t];
age_color_idx = [zeros(size(y_pred_wt)); ones(size(y_pred_a53t))];
min_age_diff = min(diff(unique(y_test)));
age_x_vals = ((y_test + 1)/min_age_diff) * 2;

wt_pred_imp = mean(vertcat(wt_age.predImp));
a53t_pred_imp = mean(vertcat(a53t_age.predImp));
[top_wt_vals, top_wt_idx] = maxk(wt_pred_imp,N_PRED);
wt_preds = age_feature_names(top_wt_idx); wt_preds = ["nTRS","GEC","sTEA","DEC","nRFT"];
[top_a53t_vals, top_a53t_idx] = maxk(a53t_pred_imp,N_PRED);
a53t_preds = age_feature_names(top_a53t_idx); a53t_preds = ["nSFR","PAF","nCCD","sPDE","nCFS"];

%% PLOT
f3b = figure('Color','w','Position',[1200 100 595 595],'Units','points');

%Heatmap subplot
s_hmap = subplot('Position',[x_loc,...
                    y_loc,...
                    size(heatmap_log_ratio,1)*boxsize,...
                    size(heatmap_log_ratio,2)*boxsize]);
hmap = imagesc(heatmap_log_ratio');
hold on
%Mark plotted features in heatmap
arrayfun(@(x) rectangle('Position',[0.5 x-0.5 5 1],'FaceColor','none','LineWidth',1.5),example_idx);
s_hmap.TickLength = [0 0];
s_hmap.FontSize = fontsz;
s_hmap.YTick = 1:length(labels);
s_hmap.YTickLabel = labels;
s_hmap.XTick = 1:5;
s_hmap.XLabel.String = "Week";
xtickangle(s_hmap,0)

colormap(heatmap_colormap)
set(gca,'CLim',[log(1/3) log(3)])
cb = colorbar('Position',[x_loc s_hmap.Position(2)+s_hmap.Position(4)+diff_colorbar, s_hmap.Position(3), cb_height],...
    'Location','northoutside');
cb.Ticks = [-1.05 0 1.05];
cb.TickLabels = ["<0.3", "1", ">3"];
cb.Ruler.TickLabelGapOffset = -1;
cb.Ruler.TickLabelRotation = 0;
cb.Title.String = "\bf WT/A53T";
cb.Title.Position = [21 14 0];

%Add feature importance scatter
[Y,X] = ndgrid(1:size(feature_importances,1),1:size(feature_importances,2));
scatter(X(:),Y(:),feature_importances(:)*10,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.8)

% Accuracies
s_fi = subplot('Position',[x_loc+size(heatmap_log_ratio,1)*boxsize+diff_heatmaps,...
    y_loc,...
    size(heatmap_log_ratio,1)*boxsize,...
    size(heatmap_log_ratio,2)*boxsize]);
b = barh(s_fi,1:length(accuracy_values),accuracy_values(end:-1:1),1,'k','EdgeColor','w','BarWidth',0.01);
hold on
% er = errorbar(b.YData,b.XData,zeros(1,length(accuracy_sd)),accuracy_sd(end:-1:1),'horizontal','Color',[0.5 0.5 0.5],...
%     'linestyle','none','CapSize',1);
title('RF','FontWeight','normal','FontSize',fontsz)
axis tight
yticks([])
xtickangle(s_fi,0)
xlabel("Accuracy")
set(gca,'FontSize',7)
box off
axis tight
set(gca,'XLim',[0.5 1])
set(gca,'XTick',[0.5, 0.75, 1.0])



%Statistical significance
s_ss = subplot('Position',[s_fi.Position(1)+s_fi.Position(3)+0.015,...
    y_loc,...
    0.05,...
    s_fi.Position(4)]);
plot(ones(1,sum(num_stars>0)),length(labels)-find(num_stars)+1,'k*','MarkerSize',3)
hold on
plot(ones(1,sum(num_stars>1))+0.1,length(labels)-find(num_stars>1)+1,'k*','MarkerSize',3)
plot(ones(1,sum(num_stars>2))+0.2,length(labels)-find(num_stars>2)+1,'k*','MarkerSize',3)
ylim([0.5 length(labels)+0.5])
xlim([1 1.5])
t_MLM = title('LMM','FontWeight','normal','FontSize',fontsz,'HorizontalAlignment','right');
t_MLM.Position(1) = 1.28;
axis off
annotation(f3b,'rectangle',...
    [0.315243697478992 0.300840336134454 0.0410168067226891 0.62],...
    'LineStyle',':','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.1);

% Draw gap between feature groups
white_gap = annotation(f3b,'rectangle',...
    [0.06 0.601 0.3 0.017],...
    'LineStyle','none',...
    'FaceColor',[1 1 1]);
feature_sep1 = yline(s_hmap,TAKE_N_BEST+0.5,'LineWidth',1);
feature_sep2 = yline(s_hmap,TAKE_N_BEST+1.5,'LineWidth',1);


%Example timeline plots
c = heatmap_colormap([1,end],:);
for p = 1:length(example_idx)
    exp_plots{p} = subplot('Position',...
        [white_gap.Position(1)+white_gap.Position(3)+0.09,...
        y_loc+(p-1)*(y_plot+0.078),...
        x_plot,...
        y_plot]);
    for i = size(timeline_data,1):-1:1
        x = tps+1+jitter(i);
        y = timeline_data(i,:,example_idx(p));
        xx = linspace(min(x),max(x),100); yy = pchip(x,y,xx);
        plot(xx,yy,'Color',c(i,:))
        hold on
        errorbar(tps+1+jitter(i),timeline_data(i,:,example_idx(p)),timeline_sd(i,:,example_idx(p)),...
            'LineWidth',1,'Color',c(i,:),'CapSize',0,'LineStyle','none','Marker','o','MarkerSize',2,...
            'MarkerFaceColor',c(i,:),'HandleVisibility','off');
        set(gca,'FontSize',fontsz)
        marg = get(gca,'ylabel');
        set(marg,'Margin',3)
    end
    box off
    if p == 1
        xlabel('Week')
    elseif p == 4
        leg = legend(["WT","A53T"],'box','off');
        leg.ItemTokenSize = [5 5];
        leg.Position = [0.5, 0.92, 0.065, 0.0445];
        leg.NumColumns = 2;
    end
    ylabel(example_labels(p))
    xlim([3 max(xlim)])
    xticks(7:7:35)
    xticklabels(1:5)
    yt = get(gca,'YTick');
    yticks(linspace(min(yt),max(yt),2)) %Change here to insert more ticks
    yl = get(gca,'YLim');
    ylim([yl(1)-max(abs(yl))*0.1 yl(2)])
end

umap_sc = subplot('Position',...
    [exp_plots{1}.Position(1)+exp_plots{1}.Position(3)+0.1,...
    y_loc+2*(y_sz+0.1),...
    x_sz,...
    y_sz]);

RecordingGroup.plot_cluster_outlines(sc_reduction, sc_true_idx, umap_sc);

title("Single-cell features" + newline + "Purity: " + sprintf('%.2f', sc_clust_coef),'FontWeight','normal')
set(gca,'FontSize',fontsz)

ump_nw = subplot('Position',...
    [exp_plots{1}.Position(1)+exp_plots{1}.Position(3)+0.1,...
    y_loc+(y_sz+0.1),...
    x_sz,...
    y_sz]);

RecordingGroup.plot_cluster_outlines(nw_reduction, nw_true_idx, ump_nw);

title("Network features" + newline + "Purity: " + sprintf('%.2f', nw_clust_coef),'FontWeight','normal')
set(gca,'FontSize',fontsz)

umap_comb = subplot('Position',...
    [exp_plots{1}.Position(1)+exp_plots{1}.Position(3)+0.1,...
    y_loc,...
    x_sz,...
    y_sz]);

RecordingGroup.plot_cluster_outlines(comb_reduction, comb_true_idx, umap_comb);

title("Combined features" + newline + "Purity: " + sprintf('%.2f', comb_clust_coef),'FontWeight','normal')
set(gca,'FontSize',fontsz)

% Heatmap groups x time
s_hgg = subplot('Position',[0.2 0.05 0.1 y_sz]);
imagesc(s_hgg,test_acc([5:8,1:4,9],:))
hold on
yline(4.5,'LineWidth',1,'Color','k')
yline(8.5,'LineWidth',1,'Color','k')
xline(5.5,'LineWidth',1,'Color','k')
xticks(1:size(test_acc,2))
xticklabels([string(1:5),"all"])
xtickangle(s_hgg,0)
xlabel('Week')
yticks(1:size(test_acc,1))
ytl = feature_groups([5:8,1:4,9]);
ytl([3,5]) = "Time series";
yticklabels(ytl)
colormap(s_hgg,greens)
cb = colorbar('Limits',[0.5 1],'Ticks',[0.5:0.1:1]);
set(s_hgg,'FontSize',fontsz)
t = title('Classification accuracy','FontWeight','normal');t.Position(2) = 0;
t.FontSize = fontsz;
s_hgg.Position(3) = 0.1;

% Age prediction
% s_ap = subplot('Position',[s_hgg.Position(1)+s_hgg.Position(3)+0.15 s_hgg.Position(2) x_sz y_sz]);
s_ap = subplot('Position',[0.44 s_hgg.Position(2) x_sz y_sz]);
bc = boxchart(s_ap, age_x_vals, y_pred,'GroupByColor',age_color_idx,'MarkerSize',1,'BoxWidth',0.5);
hold on
p = plot(s_ap, [0 max(age_x_vals)*1.1],[0 max(y_test)*1.1],'k--'); p.Color(4) = 0.3;
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xticks((1:length(tps)) * 2)
xticklabels(1:length(tps))
xlabel('Culture age [weeks]')
yticks(tps)
yticklabels(1:length(tps))
ylabel('Predicted age [weeks]')
set(s_ap,'FontSize',fontsz)
leg = legend("WT", "A53T",'Box','off','Location','northwest'); leg.ItemTokenSize = [5 5];

% s_ap_pi = subplot('Position',[s_ap.Position(1)+s_ap.Position(3)+0.1 s_ap.Position(2) x_sz y_sz]);
s_ap_pi = subplot('Position',[0.73 s_ap.Position(2) x_sz y_sz]);
pi_b = bar([top_wt_vals; top_a53t_vals],'FaceColor','flat','BarWidth',0.8);
xticklabels(["WT", "A53T"])
ylabel('Feature importance')
arrayfun(@(x) set(x,'CData', flipud(c)),pi_b)
box off
xtips = arrayfun(@(x) x.XEndPoints,pi_b,'un',0); [xtips, x_order] = sort([xtips{:}],'ascend');

ytips = arrayfun(@(x) x.YEndPoints,pi_b,'un',0); ytips = [ytips{:}];
ytips = ytips(x_order) +0.1;
set(s_ap_pi,'FontSize',fontsz)
text(xtips,ytips,[wt_preds, a53t_preds],'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Rotation',90, 'FontSize',fontsz - 0.5)
setAllFontSizes(f3b,fontsz)

%%
exportgraphics(f3b,'/home/phornauer/Git/DeePhys/Plots/Figure3b.tif','Resolution',300)