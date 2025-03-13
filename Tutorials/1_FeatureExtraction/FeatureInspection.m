% After the feature extraction we often want to look at the output and inspect
% the extracted features to identify any problems or simply for nice plots.

% For this example, we will pick a random culture and show some common
% functions to inspect the output. In practice, you probably want to do
% this a bit more systematically, especially if you are working with e.g. a
% new cell line. In particular the burst parameters are often tricky to
% infer automatically, so some cultures might need some individual
% parameter tuning. Since the cultures we are investigating are
% unforunately not bursty, I will only provide some code snippets that
% will, however, not run.


% Specify the DeePhys root path
addpath(genpath("/home/phornauer/Git/DeePhys")) % CHANGE

% Or if your working directory is 2_FeatureExtraction we can also use a
% relative path
addpath(genpath("../../Functions"))
addpath(genpath("../../Classes"))
addpath(genpath("../../Toolboxes"))

%% Load data
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/"; % Where you saved the data to
file_path = "iNeurons_dataset/network/240610/T002513/Network/well015/sorter_output/MEArecording.mat"; % Select another object by changing the well ID (e.g. well002)
full_path = fullfile(root_path,file_path);

load(full_path)

%% The MEArecording obj
% First, take a look at the Properties of our 'obj'. It should (hopefully)
% be pretty much self explanatory.

obj

% Metadata: Inferred from the lookup.xlsx file you provided during feature
%   extraction
% RecordingInfo: Info inferred from the sorting output (duration etc)
% Parameters: Provided at the feature extraction, otherwise default values
%   that you can also find under obj.returnDefaultParams()
% Spikes: Spike times/units from those units that passed the QC
% Bursts: Structure with information about the automated burst detection.
%   Preserves the raw burst detection output and the curated one (not
%   available for the tutorial dataset)
% Units: Units that passed the QC stored as a separate object. Contains all
%   information about the unit.
% NetworkFeatures: Contains network feature groups that in return contain
%   tables with the respective extracted feature values.
% UnitFeatures: Similar to NetworkFeatures, but contains the aggregate of
%   all unit feature values across the recording.
% ClusteredFeatures: Similar to UnitFeatures, but is calculated after
%   running a clustering algorithm. Instead of one value for the entire
%   culture, this will give one value for each cluster (see
%   SingleCellPhenotype tutorial).
% Connectivity: Contains graphs inferred from the connectivity algorithms
%   and the CCGs (if selected during analysis).

%% Retrieve relevant data from the obj
% Some convenience functions :

% To retrieve the 'original' spike sorting spike times
[spike_times, spike_units, hasSpikes] = getSpikeTimes(obj);

% Or the 'original' template matrix
template_matrix = getTemplateMatrix(obj);

% Retrieve a table of all single-cell/unit features for all units
unit_table = obj.getUnitFeatures(obj,"all");

%% Example use case: Looking for local connectivity
% First, we look for active neurons. This is particulary important when the
% recordings were short, so that we have enough activity for proper CCGs.
[fr, fr_idx] = maxk(unit_table.FiringRate, 10);

% Then we look for units that were found have a significant connection
% using the CCG method.
exc_con = obj.Connectivity.CCG.ExcitatoryConnection;

% Then we look for active units in the connection list.
unit_overlap = ismember(obj.Connectivity.CCG.ExcitatoryConnection, fr_idx);

% We choose the one where both units are among the most active.
candidate_con = find(sum(unit_overlap,2) == 2);

% Let's check the actual CCG
figure('Color','w')
obj.PlotCCG(exc_con(candidate_con,1),exc_con(candidate_con,2)); 
% Great, we see a nice peak! Now we can additionally check where the units
% are located to exclude e.g. spike sorting errors.

unit_ids = [exc_con(candidate_con,1),exc_con(candidate_con,2)];
wf_cutout = 1:60; %Samples from the waveform to plot
scale_factor = 1.9; %Needs to be adjusted depending on the distance of the electrodes
figure('Color','w')
obj.PlotUnitTemplate(unit_ids, scale_factor, wf_cutout);

% Let's overlay the plot with another nearby unit
min_dist = 25; %Minimum distance from AIS to AIS 
max_dist = 150; %Maximum distance from AIS to AIS 
close_units = obj.findCloseUnits(exc_con(candidate_con,1), min_dist, max_dist);

unit_ids = [exc_con(candidate_con,1),exc_con(candidate_con,2), close_units];
figure('Color','w')
obj.PlotUnitTemplate(unit_ids, scale_factor, wf_cutout);

%% Some additional functions regarding connectivity
% Obtain the empirical STTC matrix
empirical_sttc = obj.empiricalSTTCmatrix();

% Plot connectivity matrix
method = "STTC";
matrix = "wu"; % Weighted undirected connectivity matrix
figure('Color','w')
obj.PlotCommunity(method,matrix);

% Cluster activity (needs clustering from SingleCellPhenotype tutorial)
time_cutout = [0, 100];
% obj.PlotUnitClusterActivity(time_cutout); %This would return an error without the clustering results

%% Working with bursts
% For demonstration purposes only, we will go through the burst detection
% workflow here. Since there were no bursts in the cultures, we do not
% expect to find any meaningful data here.
obj.detectBursts();

% As expected, we didn't find any bursts, let's check if that's correct.
obj.PlotBurstCheck("hist",cutout=[150,200],binning=0.1)
% The spike in coactivity might be burst-like, let's look at the scatter
% plot.
obj.PlotBurstCheck("scatter",cutout=[150,200])
% Seems like the burst detection did a good job in not labeling anything as
% a burst.

%% Automatic burst detection 
% In short, the algorithm will infer the best ISI_N and N values based on
% the trough in the ISI distribution (see original paper) and the dunns
% coefficient. This is somewhat heuristic, so if you are not satisfied with
% the results, you can also specify those two values directly. 

% Checking the parameter combinations suggested by the algorithm is mostlz
% still a good starting point:
obj.Bursts
% obj.Bursts.N and obj.Bursts.ISI_N are the parameters that were tested,
% obj.Bursts.dunns gives the corresponding dunn's coefficient (higher is
% better). 

%% 'Manual'burst detection 
% Let's try different parameters to see what happens (and to showcase some 
% additional functionality):

obj.detectBursts(N=7, ISI_N=0.01);

% Only small differences from the inferred parameters already changes the
% result to several detected bursts. So if you do not expect bursts it is best
% not to run the burst detection in the first place, as the probability of
% a few false positives is pretty high.

% The results from the original burst detection are stored in
obj.Bursts.raw
% which we can also use for the plot function
obj.PlotBurstCheck("hist",cutout=[50,55],binning=0.1,input="raw")

% Those bursts are then pruned and merged. If you feel like the bursts are
% overly split or merged, you can adjust the input parameter 'merge_factor' 
% to the detectBursts() function.

%% Saving changes
% If you made changes that you want to persist after you close your current
% MATLAB session:

% The overwrite flag needs to be set (should be the default)
obj.Parameters.Save.Overwrite = true;
obj.saveObject();
