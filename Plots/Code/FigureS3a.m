addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'SCR*','230113','*','w*','sorted'}; 
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {{'Source','Taylor'},{'Mutation',"WT"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
%,{'Mutation','WT'}
rg_params.Selection.Exclusion = {{'DIV',12},{'Treatment',"ASO","LNA"}};%,{'Treatment',"LNA"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt = RecordingGroup(rec_array, rg_params);

wt_units = arrayfun(@(x) length(x.Units), wt.Recordings);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {{'Source','Taylor'},{'Mutation',"SNCA"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
%,{'Mutation','WT'}
rg_params.Selection.Exclusion = {{'DIV',12},{'Treatment',"ASO","LNA"}};%,{'Treatment',"LNA"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
snca = RecordingGroup(rec_array, rg_params);

snca_units = arrayfun(@(x) length(x.Units), snca.Recordings);

%% Scatter plot
reds = othercolor('Reds9',9);
red = reds(7,:);
blues = othercolor('Blues9',9);
blue = blues(7,:);
fontsz = 6.5;

time_cutout = [60 120];
fs3a = figure('Color','w','Position',[1200 100 595 595],'Units','points');
tl = tiledlayout(7,4,'TileSpacing','compact','Padding','tight');
nt(1) = nexttile(tl,1);
wt.Recordings(4).PlotNetworkScatter(time_cutout,blue)
nt(1).YTick = [100, 200, 300];
nt(1).XTick = [60, 120];
xticklabels([0,60])
ylabel('Units')
% nt(1).XLabel = [];%title('WT','Color',blue)

nt(2) = nexttile(tl,2);
snca.Recordings(3).PlotNetworkScatter(time_cutout,red)
nt(2).YTick = [100, 200];
nt(2).XTick = [60, 120];
xticklabels([0,60])
ylabel('Units')
%title('A53T','Color',red)
setAllFontSizes(fs3a, fontsz)
exportgraphics(fs3a,'/home/phornauer/Git/DeePhys/Plots/FigureS3a.tif','Resolution',300)