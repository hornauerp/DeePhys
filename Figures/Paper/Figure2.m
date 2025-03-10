%% Batch 3 FCDI
rg = RecordingGroup(fcdi_array);

%% Select most active recording and units
n_unit_spikes = arrayfun(@(x) length([x.SpikeTimes]),rg.Units);
[~, max_unit_idx] = maxk(n_unit_spikes,10);

n_nw_spikes = arrayfun(@(x) length([x.Spikes.Times]),rg.Recordings);
[~, max_nw_idx] = max(n_nw_spikes);

%% Bin spiketimes within selected range
N_seconds = 90;
binned_mat = nan(length(max_unit_idx),N_seconds * 10);

for i = 1:length(max_unit_idx)
    binned_mat(i,:) = histcounts(rg.Units(max_unit_idx(i)).SpikeTimes,0:0.1:N_seconds);
end

%% Generate distribution
binned_nw = histcounts(rg.Recordings(max_nw_idx).Spikes.Times,0:0.1:N_seconds);
nw_dist = smoothdata(histcounts(binned_nw,0:2:max(binned_nw)),'SmoothingFactor',0.05);
outlier_idx = round(length(nw_dist) * 0.7);
nw_clean = nw_dist(1:outlier_idx);
nw_outlier = nw_dist(length(nw_clean) + 1 : end);

binned_nw = histcounts(rg.Recordings(max_nw_idx).Spikes.Times,0:0.1:900);
[acf, lags] = autocorr(binned_nw,'NumLags',3000);

norm_nw = normalize(binned_nw(1:10:length(acf)),'range',[0 1]);
norm_acf = normalize(acf(1:10:end),'range',[0 1]);

%%
norm_nw_activity = binned_nw/max(binned_nw);
NFFT = length(norm_nw_activity);
F = (0 : 1/NFFT : 1/2-1/NFFT)*(1/0.1);
TEMP = fft(norm_nw_activity,NFFT);
TEMP(1) = 0;
freq_domain = abs(TEMP(1:NFFT/2));

%% Allocate colors
reds = othercolor('Reds9',9);
red = reds(7,:);
greens = othercolor('Greens9',9);
green = greens(7,:);
blues = othercolor('Blues9',9);
blue = blues(7,:);
greys = othercolor('Greys9',30);
greys = greys(16:25,:);

fontsz = 8;

%%
f2 = figure('Color','w','Position',[2478,543,224,475]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','tight');
ax(1) = nexttile(1,[1 2]);
hold on
arrayfun(@(x) plot(ax(1), 1:size(binned_mat,2),binned_mat(x,:) + max(max(binned_mat)) * (x -1),'Color',greys(x,:)),1:length(max_unit_idx))
xticks(0:300:900);xticklabels(0:30:90);xlim([0 900])
xlabel('Time [s]')
% yticks(max(max(binned_mat)) +2 : max(max(binned_mat)) : max(max(binned_mat)) * length(max_unit_idx) + 2)
ylabel('Units')
ax(1).YAxis.Visible = 'off'; ax(1).TickLength = [0 0];
ax(1).YRuler.TickLabelGapOffset = -9;
set(findall(gca, 'type', 'text'), 'visible', 'on')

ax(2) = nexttile(4);
hold on
plot(ax(2), 1:length(nw_clean), nw_clean, 'k')
plot(ax(2), length(nw_clean) + 1:length(nw_dist), nw_outlier,'Color',red)

xline(ax(2), mean(nw_dist),':','LineWidth',0.1,'Label','Mean','LabelHorizontalAlignment','center','FontSize',fontsz)
xline(ax(2), outlier_idx, ':','Color',red,'LineWidth',0.1,'Label','Outlier','LabelHorizontalAlignment','center','FontSize',fontsz)
set(gca,'XScale','log')
xlabel(ax(2),'     Coactivity\newline [spikes/100ms]')
ylabel('Count [bins]')
title('Distribution','FontWeight','normal')

% ax = subplot(2,2,4);
% plot(ax, binned_nw(1:10:end),'k','LineWidth',1)
% hold on
% plot(ax, smoothdata(binned_nw(1:10:end),'SmoothingFactor',0.5),'Color',green,'LineWidth',2,'LineStyle',':')
% yticks([])
% xticks([])
% xlabel('Time')
% box off
% title('Fluctuation analysis','FontSize',7)

ax(3) = nexttile(3);
plot(ax(3), norm_nw,'Color',greys(end,:))
hold on
plot(ax(3), norm_acf,'Color',green,'LineWidth',2,'LineStyle',':')
yticks([0 1])
xlim([0 300])
xticks([0 300])
xlabel('Time [s]')
ylabel('Normalized data')

ylim([0 1.2])
box off
title('Autocorrelation','FontWeight','normal')
leg = legend(["Activity","ACF"],'Box','off','Location','north'); leg.Position(1:2) = leg.Position(1:2) + [0.04 0.02];
leg.ItemTokenSize = [5 5];


% set(findall(gcf, 'type', 'text'), 'FontSize', fontsz)
setAllFontSizes(f2, fontsz)

%%
exportgraphics(f2,'/home/phornauer/Git/DeePhys/Plots/Figure2.tif','Resolution',300)