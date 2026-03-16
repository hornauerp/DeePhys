% CellTypeClassification.m
%
% Classify units into excitatory and inhibitory cell types using the
% @CellTypeClassifier pipeline. This leverages pharmacological perturbation
% data (e.g. DREADD/chemogenetic activation) to identify inhibitory
% candidates via bootstrap firing-rate tests, then trains a supervised UMAP
% classifier using waveform and ACG features.
%
% The pipeline:
%   1. Load recordings and form a RecordingGroup with culture grouping
%   2. Initialise CellTypeClassifier with parameters
%   3. Bootstrap firing-rate test to identify responsive (inhibitory) units
%   4. Generate training labels via UMAP + isolation forest
%   5. Classify all units with supervised UMAP
%   6. Optionally apply to external datasets
%
% This replaces the old manual workflow in clf_transfer.m / clf_all_batches.m
% which used standalone functions (identify_interneurons, bootstrap_fr_response,
% etc.) and required manual stitching of the pipeline steps.
%
% REQUIREMENTS:
%   - DeePhys on path
%   - run_umap toolbox (MATLAB File Exchange) for supervised UMAP

%% 0 - Setup
addpath(genpath("/path/to/DeePhys")) % CHANGE

%% 1 - Load recordings
root_path = "/path/to/your/data/"; % CHANGE
path_logic = {'T00*', 'Network', 'w*', 'sorter_output'};
path_list = generate_sorting_path_list(root_path, path_logic);

rm_con = false; % keep connectivity for this workflow
mearec_array = recording_array_from_single_files(path_list, rm_con);

%% 2 - Form RecordingGroup
% Filter to recordings of interest. Example: cortical cultures with
% specific AAV constructs for chemogenetic experiments.
rg_params.Selection.Inclusion = {{'Region', 'Cortex'}, {'AAV', 0, 128, 129}};
rg_params.Selection.Exclusion = {{'Concentration', inf}};
rg = RecordingGroup(mearec_array, rg_params);

% Group into cultures by dose-response series (one culture = one chip
% across concentrations).
rg.Cultures = rg.groupCultures('RecordingDate');

%% 3 - Initialise CellTypeClassifier
parameters = struct();

% Bootstrap firing-rate test
parameters.Bootstrap.PreCutout  = [0, 1200];    % baseline window (s)
parameters.Bootstrap.PostCutout = [6000, 7200]; % post-treatment window (s)
parameters.Bootstrap.BinSize    = 20;           % bin width (s)
parameters.Bootstrap.NIter      = 1000;
parameters.Bootstrap.Alpha      = 1e-10;

% UMAP embedding
parameters.UMAP.UnitFeatures    = ["FullACG", "ReferenceWaveform"];
parameters.UMAP.NormalizationVar = "ChipID";
parameters.UMAP.GroupingVar     = "Concentration";
parameters.UMAP.GroupingValues  = 0;            % baseline only
parameters.UMAP.NNeighbors      = 100;
parameters.UMAP.NDims           = 10;

% Outlier detection (isolation forest on inhibitory candidates)
parameters.OutlierDetection.ContaminationFraction = 0.5;
parameters.OutlierDetection.NObsPerLearner        = 50;
parameters.OutlierDetection.DistancePercentile    = 80;

ctc = CellTypeClassifier(rg, parameters);

%% 4 - Identify inhibitory candidates
% For cultures with a specific treatment (e.g. AAV=128 = DREADD-inhibitory),
% test each unit for significant firing-rate change between pre and post.
metadata_filter = {'AAV', 128};
ctc.identifyResponsiveUnits(metadata_filter);

fprintf('Responsive units: %d/%d\n', ...
    sum(ctc.ResponsiveUnitIdx), length(ctc.ResponsiveUnitIdx))

%% 5 - Generate training labels
% Runs unsupervised UMAP on all units, then uses isolation forest to
% refine the inhibitory cluster and select counterexamples (excitatory).
ctc.generateTrainLabels();

% Training labels:
%   ctc.TrainLabels.sorted_y_train — 1=excitatory, 2=inhibitory
%   ctc.TrainLabels.sorted_train_ids — indices into ctc.UnitList

%% 6 - Classify all units
% Supervised UMAP projects train + test units, propagates labels.
ctc.classifyUnits();

% Results:
%   ctc.UnitLabels — 1=excitatory, 2=inhibitory, NaN=unclassified
n_exc = sum(ctc.UnitLabels == 1);
n_inh = sum(ctc.UnitLabels == 2);
fprintf('Classification: %d excitatory, %d inhibitory\n', n_exc, n_inh)

%% 7 - Inspect classification quality
% Check that inhibitory units have distinct ACG and waveform profiles
figure('Color', 'w');
subplot(1, 2, 1)
title('Supervised UMAP')
scatter(ctc.Reduction.Test(:,1), ctc.Reduction.Test(:,2), 5, ctc.UnitLabels, 'filled')
colorbar

% Per-chip I/E ratio sanity check
chips = unique([cellfun(@(x) x(1).Metadata.ChipID, rg.Cultures, 'un', false)]);
for c = 1:min(5, length(chips))
    fprintf('Chip %s: ', chips{c})
    % count E/I per chip
end

%% 8 - Apply to external data (optional)
% If you have waveforms and ACGs from a different recording system:
%
% external_waveforms — (N_units x N_samples) matrix
% external_acgs      — (N_units x N_bins) matrix
%
% result = ctc.classifyExternalUnits(external_waveforms, external_acgs);
% external_labels = result.Labels; % 1=exc, 2=inh

%% 9 - Alternative: metadata-based labels
% If you dont have pharmacological data but have metadata labels (e.g.
% from manual curation or a separate experiment), you can set train labels
% directly:
%
% ctc.TrainLabels.sorted_y_train = your_labels;
% ctc.TrainLabels.sorted_train_ids = your_unit_indices;
% ctc.classifyUnits();

%% 10 - Save results
% The CellTypeClassifier object contains all intermediate results.
% Save it for downstream EIAnalysis:
% save('ctc_results.mat', 'ctc');
