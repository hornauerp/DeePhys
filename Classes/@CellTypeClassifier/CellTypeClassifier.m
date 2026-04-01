classdef CellTypeClassifier < handle
% CELLTYPECLASSIFIER  Supervised cell-type classification pipeline for MEA experiments.
%
% Identifies inhibitory (interneuron) vs excitatory neurons using a bootstrap
% firing-rate response test followed by supervised UMAP classification.
%
% USAGE (new API — v2):
%   ctc = CellTypeClassifier(featureStore, unitDataArray, params);
%   ctc.identifyResponsiveUnits();    % bootstrap pre/post firing rate comparison
%   ctc.generateTrainLabels();        % UMAP embedding + isolation forest label generation
%   ctc.classifyUnits();              % supervised UMAP projection and classification
%   labels = ctc.UnitLabels;          % 1 = excitatory, 2 = inhibitory, NaN = unclassified
%
% MIGRATION from old API:
%   ctc = CellTypeClassifier.fromLegacyGroup(rg, params);

    properties
<<<<<<< HEAD
        FeatureStore            FeatureStore    % Features + metadata for all units
        UnitDataArray           UnitData        % (1 x N) raw unit data (spike times etc.)
        Parameters              struct          % Merged from returnDefaultParams + user overrides
        ResponsiveUnitIdx       logical         % (1 x N) units with a significant firing rate response
        ResponsiveUnitDirection string          % (1 x N) "none" | "increase" | "decrease" per unit
        TrainLabels             struct          % .sorted_train_ids, .sorted_y_train, etc.
        UnitLabels              double          % (1 x N): 1=excitatory, 2=inhibitory, NaN=unclassified
        UnitConfidence          double          % (1 x N): kNN confidence in [0,1]; 1.0 for training units
        UMAP                                   % Trained UMAP model from run_umap
        Reduction               struct          % Struct with Unsupervised/Train/Test/External embeddings
        HarmonizedWaveforms     double          % (N_samples × N_units) processed waveforms
        HarmonizedACGs          double          % (N_bins × N_units) ACGs
        HarmonizedSR            double          % Waveform sampling rate after harmonization
        CachedExtraction        struct          % Cached waveform/ACG extraction
    end

    % Kept for backward compatibility (populated by fromLegacyGroup)
    properties (Access = private)
        LegacyRecordingGroup    % RecordingGroup, only set by fromLegacyGroup
=======
        RecordingGroup          % RecordingGroup with .Cultures already set
        Parameters              % struct — merged from returnDefaultParams + user overrides
        UnitList                % (1 x N) Unit array — all units across all cultures
        ResponsiveUnitIdx       % (1 x N) logical — units identified as responsive (inhibitory candidates) by bootstrap
        ResponsiveStrength      % (1 x N) double — continuous response strength per unit (Spearman rho or effect size)
        TrainLabels             % struct: .sorted_train_ids, .sorted_y_train, .umap_train_idx, .umap_test_idx
        UnitLabels              % (1 x N) double — 1=excitatory, 2=inhibitory, NaN=unclassified
        UMAP                    % trained UMAP model returned by run_umap
        Reduction = struct('Unsupervised', [], 'Train', [], 'Test', [], 'External', [], 'Extras', [])
        HarmonizedWaveforms     % (N_samples × N_units) processed waveforms (upsampled, aligned, trimmed)
        HarmonizedACGs          % (N_bins × N_units)   ACGs at Harmonization params
        HarmonizedSR            % scalar — waveform sampling rate after harmonization
        NormalizationParams     % struct caching normalization from generateTrainLabels for carryover to classifyUnits/applyTo
                                %   .mu_global    — (1 x F) mean from global z-score
                                %   .sigma_global — (1 x F) std from global z-score
                                %   .nan_cols     — (1 x F) logical mask of NaN columns removed
                                %   .scale        — (1 x F') max(abs) scaling vector (after NaN col removal)
        Confidence              % (1 x N) double — per-unit confidence score (0-1, 1=most confident)
        CachedExtraction        % struct caching extractUnitWaveformsAndACGs output (see getOrExtract)
                                % struct storing UMAP embeddings:
                                %   .Unsupervised — (N x D) from generateTrainLabels (all units, unsupervised)
                                %   .Train        — (N_train x D) from classifyUnits (supervised training embedding)
                                %   .Test         — (N_test  x D) from classifyUnits (supervised test projection)
                                %   .External     — (N_ext  x D) from classifyExternalUnits
>>>>>>> claude/dreamy-varahamihira
    end

    methods

        function ctc = CellTypeClassifier(feature_store, unit_data, parameters)
        % CELLTYPECLASSIFIER  Construct from FeatureStore + UnitData array.
        %
        %   ctc = CellTypeClassifier(featureStore, unitDataArray)
        %   ctc = CellTypeClassifier(featureStore, unitDataArray, params)
        %
        % unit_data provides spike times for the bootstrap test; feature_store
        % provides waveform/ACG features and metadata for UMAP.
            arguments
                feature_store   FeatureStore = FeatureStore()
                unit_data       UnitData     = UnitData.empty()
                parameters      struct       = struct()
            end
            if isempty(feature_store.UnitTable)
                return
            end
            assert(numel(unit_data) == height(feature_store.UnitTable), ...
                ['CellTypeClassifier: UnitDataArray length (%d) must match ' ...
                 'UnitTable height (%d). Build both from the same FeatureStore ' ...
                 'assembly call.'], numel(unit_data), height(feature_store.UnitTable));
            ctc.FeatureStore  = feature_store;
            ctc.UnitDataArray = unit_data;
            ctc.Parameters    = parseStructParameters(ctc.returnDefaultParams(), parameters);
        end

        function [wf, acg, sr] = getOrExtract(ctc, unit_data)
        % GETOREXTRACT  Three-level cache: in-memory → disk → compute.
        %
        % ACG source strategy (ACGSource parameter):
        %   "FullACG" — prefer Parent_ACG* columns from FeatureStore (pre-computed
        %               from full parent spike train). Falls back to segment ACG if
        %               Parent_ACG* columns are not present in the FeatureStore.
        %   "ACG"     — always compute from UnitData.SpikeTimes (segment-level).
            ph    = ctc.Parameters.Harmonization;
            N     = numel(unit_data);
            cache = ctc.CachedExtraction;

            % ── Level 1: in-memory cache ──────────────────────────────────────
            if ~isempty(cache) ...
                    && cache.N == N ...
                    && abs(cache.ACGBinSize - ph.ACGBinSize) < 1e-12 ...
                    && abs(cache.ACGLag - ph.ACGLag) < 1e-12 ...
                    && cache.ACGSource == ph.ACGSource
                wf  = cache.wf;
                acg = cache.acg;
                sr  = cache.sr;
                return
            end

            % ── Level 2: disk cache ───────────────────────────────────────────
            % Unit identity encoded as a compact hash of sorted IDs (avoids
            % hashing a large cell array for datasets with many units).
            uid_hash   = paramHash(strjoin(sort(string({unit_data.UnitID})), '|'));
            key_struct = struct('ACGBinSize', ph.ACGBinSize, 'ACGLag', ph.ACGLag, ...
                'ACGSource', ph.ACGSource, 'N', N, 'UnitIDHash', uid_hash);
            cache_file = CacheManager.instance().get(paramHash(key_struct));

            if exist(cache_file, 'file')
                s   = load(cache_file, 'wf', 'acg', 'sr');
                wf  = s.wf;
                acg = s.acg;
                sr  = s.sr;
            else
                % ── Level 3: compute ──────────────────────────────────────────
                % Waveform: always from UnitData
                wf = double([unit_data.ReferenceWaveform]);
                if isempty(wf)
                    wf = zeros(0, N);
                    sr = ph.WaveformTargetSamplingRate;
                else
                    sr = unit_data(1).SamplingRate;
                end

                % ACG: prefer Parent_ACG* from FeatureStore when ACGSource == "FullACG"
                acg = ctc.resolveACG(unit_data, ph);

                % Write to disk cache
                try
                    save(cache_file, 'wf', 'acg', 'sr', '-v7.3');
                catch
                end
            end

            % ── Update in-memory cache ────────────────────────────────────────
            ctc.CachedExtraction = struct('wf', wf, 'acg', acg, 'sr', sr, ...
                'N', N, 'ACGBinSize', ph.ACGBinSize, 'ACGLag', ph.ACGLag, ...
                'ACGSource', ph.ACGSource);
        end

        function acg = resolveACG(ctc, unit_data, ph)
        % RESOLVEACG  Return ACG matrix, preferring Parent_ACG* from FeatureStore.
            if ph.ACGSource == "FullACG" && ~isempty(ctc.FeatureStore)
                ut       = ctc.FeatureStore.UnitTable;
                all_cols = string(ut.Properties.VariableNames);
                parent_acg_cols = all_cols(startsWith(all_cols, "Parent_ACG"));
                if ~isempty(parent_acg_cols)
                    % Map unit_data order to FeatureStore row order by UnitID
                    ud_ids    = string({unit_data.UnitID})';
                    fs_ids    = ut.UnitID;
                    [~, loc]  = ismember(ud_ids, fs_ids);
                    valid     = loc > 0;
                    n_bins    = numel(parent_acg_cols);
                    acg       = zeros(n_bins, numel(unit_data));
                    if any(valid)
                        acg(:, valid) = ut{loc(valid), parent_acg_cols}';
                    end
                    return
                end
                fprintf('CellTypeClassifier: Parent_ACG not in FeatureStore — using segment ACG.\n');
            end
            % Fall back to segment-level computation from SpikeTimes
            [~, acg, ~] = ctc.extractFromUnitData(unit_data);
        end

    end

    % =====================================================================
    % Static constructors and defaults
    % =====================================================================
    methods (Static)

        function ctc = fromLegacyGroup(rg, parameters)
        % FROMLEGACYGROUP  Migrate from an old RecordingGroup.
            arguments
                rg          RecordingGroup
                parameters  struct = struct()
            end
            fs  = FeatureStore.fromLegacyRecordings(rg.Recordings);
            ud  = CellTypeClassifier.unitDataFromRecordingGroup(rg);
            ctc = CellTypeClassifier(fs, ud, parameters);
            ctc.LegacyRecordingGroup = rg;  % keep for any legacy method paths
        end

        function defaultParams = returnDefaultParams()
            defaultParams.Bootstrap.PreCutout   = [0, 1200];
            defaultParams.Bootstrap.PostCutout  = [6000, 7200];
            defaultParams.Bootstrap.BinSize     = 20;
            defaultParams.Bootstrap.NIter       = 1000;
            defaultParams.Bootstrap.Alpha       = 1e-10;
            defaultParams.Bootstrap.Direction   = "increase";  % "increase" | "decrease" | "both"

            defaultParams.UMAP.NDims            = 10;
            defaultParams.UMAP.NNeighbors       = 50;
            defaultParams.UMAP.MinDist          = 0.1;
            defaultParams.UMAP.Spread           = 1;
            defaultParams.UMAP.SupervisedNDims      = 2;
            defaultParams.UMAP.SupervisedNNeighbors = 100;
            defaultParams.UMAP.UnitFeatures     = ["FullACG","ReferenceWaveform"];
            defaultParams.UMAP.NormalizationVar = "ChipID";
            defaultParams.UMAP.GroupingVar      = "Concentration";
            defaultParams.UMAP.GroupingValues   = 0;
            defaultParams.UMAP.TrainingCultureIdx = [];
            defaultParams.UMAP.TemplateDir      = tempdir;
            defaultParams.UMAP.ConfidenceK      = 15;  % k for kNN confidence scoring

            defaultParams.TrainLabels.ResponsiveClassLabel = 2;  % 2=responsive→inhibitory, 1=responsive→excitatory

            defaultParams.OutlierDetection.OutlierAlpha           = 0.01;   % chi-squared / posterior threshold
            defaultParams.OutlierDetection.MaxResponsiveComponents = 3;     % max GMM components for multimodal responsive populations
            defaultParams.OutlierDetection.DistancePercentile     = 80;
            defaultParams.OutlierDetection.CounterexampleRatio    = 2;

            defaultParams.RNGSeed     = 42;
            defaultParams.CultureKeys = ["ChipID", "PlatingDate"];

            defaultParams.Harmonization.WaveformTargetSamplingRate = 120000;
            defaultParams.Harmonization.WaveformPreTrough          = 1.0;
            defaultParams.Harmonization.WaveformPostTrough         = 1.0;
            defaultParams.Harmonization.WaveformEdgeMode           = "trim";
            defaultParams.Harmonization.ACGBinSize                 = 0.0005;
            defaultParams.Harmonization.ACGLag                     = 0.1;
            defaultParams.Harmonization.ACGSource                  = "FullACG";
        end

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Access = private)

        function [wf, acg, sr] = extractFromUnitData(ctc, unit_data)
        % Extract raw waveforms and ACGs from UnitData array.
        % Returns un-harmonized waveforms at the original sampling rate.
        % buildFeatureMatrix handles all harmonization (interpolation,
        % alignment, normalization) so we avoid double-processing.
            ph  = ctc.Parameters.Harmonization;
            N   = numel(unit_data);

            % Waveform: return raw at original sampling rate
            raw_wf = [unit_data.ReferenceWaveform];  % (n_samp x N)
            if isempty(raw_wf)
                wf = zeros(0, N);
                sr = ph.WaveformTargetSamplingRate;  % placeholder when no waveforms
            else
                wf = double(raw_wf);
                sr = unit_data(1).SamplingRate;       % original sampling rate
            end

            % ACG: compute from spike times at the specified bin size / lag
            n_bins = round(ph.ACGLag / ph.ACGBinSize) * 2 + 1;
            acg    = zeros(n_bins, N);
            for u = 1:N
                acg(:, u) = CellTypeClassifier.computeACG( ...
                    unit_data(u).SpikeTimes, ph.ACGBinSize, ph.ACGLag);
            end
        end

    end

    methods (Static)
<<<<<<< HEAD

        function binned = binnedSpikeMatrix(unit_data, time_window, bin_size)
        % BINNEDSPIKEMATRIX  (N_units × N_bins) spike count matrix from UnitData.SpikeTimes.
        %   Replaces Culture.getBinnedSpikeMat without requiring Culture objects.
        %   Public so that EIAnalyzer methods can call it directly.
            bins   = time_window(1) : bin_size : time_window(2);
            N_u    = numel(unit_data);
            binned = zeros(N_u, numel(bins) - 1);
            for u = 1:N_u
                st  = unit_data(u).SpikeTimes;
                in_win = st >= time_window(1) & st < time_window(2);
                if any(in_win)
                    binned(u, :) = histcounts(st(in_win), bins);
                end
            end
=======
        function defaultParams = returnDefaultParams()
            defaultParams.Bootstrap.GroundTruthMethod = "full_curve"; % "full_curve" | "two_window"
                % full_curve — Spearman correlation across all recordings/concentrations (recommended)
                % two_window — legacy: compare pre vs post cutout windows only
            defaultParams.Bootstrap.PreCutout   = [0, 1200];    % seconds — baseline window (two_window only)
            defaultParams.Bootstrap.PostCutout  = [6000, 7200]; % seconds — post-treatment window (two_window only)
            defaultParams.Bootstrap.MinRecordings = 4;          % minimum recordings per culture for full_curve (falls back to two_window)
            defaultParams.Bootstrap.BinSize     = 20;           % seconds per bin
            defaultParams.Bootstrap.NIter       = 1000;         % bootstrap iterations
            defaultParams.Bootstrap.Alpha       = 1e-10;        % significance level
            defaultParams.Bootstrap.MinResponsiveStrength = 0;  % minimum strength to accept (0 = no filtering)
                % full_curve: absolute Spearman rho (e.g. 0.5 requires moderate monotonicity)
                % two_window: effect size |diff|/null_std (e.g. 3.0 requires strong effect)

            defaultParams.UMAP.NDims            = 10;           % unsupervised UMAP output dimensions (label generation)
            defaultParams.UMAP.NNeighbors       = 50;           % unsupervised UMAP n_neighbors
            defaultParams.UMAP.MinDist          = 0.1;
            defaultParams.UMAP.Spread           = 1;
            defaultParams.UMAP.SupervisedNDims      = 2;        % supervised UMAP output dimensions (classification)
            defaultParams.UMAP.SupervisedNNeighbors = 100;      % supervised UMAP n_neighbors
            defaultParams.UMAP.UnitFeatures     = ["FullACG","ReferenceWaveform"]; % feature groups to use
            defaultParams.UMAP.NormalizationVar = "ChipID";     % metadata field for within-group normalisation
            defaultParams.UMAP.GroupingVar      = "Concentration"; % metadata field grouping recordings (e.g. dose)
            defaultParams.UMAP.GroupingValues   = 0;            % value(s) of GroupingVar to use (0 = baseline)
            defaultParams.UMAP.TrainingCultureIdx = [];         % culture indices for label generation UMAP (empty = all)
            defaultParams.UMAP.TargetWeight     = 0.5;          % balance between data topology (0) and label topology (1) in supervised UMAP
            defaultParams.UMAP.TemplateDir      = tempdir;      % directory for UMAP template .mat file

            defaultParams.OutlierDetection.ContaminationFraction = 0.5;
            defaultParams.OutlierDetection.NObsPerLearner        = 50;
            defaultParams.OutlierDetection.CounterexampleRatio   = 2;   % excitatory:inhibitory sampling ratio

            defaultParams.Classification.UseConfidenceThreshold  = false; % set true to mark low-confidence predictions as NaN
            defaultParams.Classification.ConfidenceThreshold     = 0.3;   % units below this confidence become unclassified (only when UseConfidenceThreshold=true)
            defaultParams.Classification.UseJoinedTransform      = false; % set true to use UMAP joined_transform (reduces false positives for divergent train/test)

            defaultParams.BayesianOptimization.MaxEvaluations     = 30;  % bayesopt iterations
            defaultParams.BayesianOptimization.NNeighborsRange     = [5, 200];
            defaultParams.BayesianOptimization.MinDistRange        = [0.01, 1.0];
            defaultParams.BayesianOptimization.SupervisedNNeighborsRange = [5, 300];
            defaultParams.BayesianOptimization.ContaminationRange  = [0.1, 0.9];
            defaultParams.BayesianOptimization.SpreadRange         = [0.5, 5.0];
            defaultParams.BayesianOptimization.TargetWeightRange   = [0.1, 0.9];
            defaultParams.BayesianOptimization.NDimsRange          = [3, 20];  % unsupervised UMAP dimensions
            defaultParams.BayesianOptimization.CounterexampleRatioRange = [1, 5];
            defaultParams.BayesianOptimization.ObjectiveMetric     = "silhouette"; % "silhouette" | "qf_dissimilarity" | "combined"
            defaultParams.BayesianOptimization.UseInterneuronPenalty = false; % set true to penalise deviation from target fraction
            defaultParams.BayesianOptimization.InterneuronTarget   = 0.19; % expected interneuron fraction (only used when UseInterneuronPenalty=true)
            defaultParams.BayesianOptimization.InterneuronWeight   = 5.0;  % penalty weight (only used when UseInterneuronPenalty=true)

            defaultParams.RNGSeed = 42;  % seed for reproducible label generation & classification

            % Harmonization parameters — shared target spec for DeePhys and external data.
            % buildFeatureMatrix and extractUnitWaveformsAndACGs both read these to
            % ensure DeePhys and external units are projected into the same feature space.
            % Change these to harmonize to a different ACG resolution or waveform rate.
            defaultParams.Harmonization.WaveformTargetSamplingRate = 120000;  % Hz — output rate after interpolation
            defaultParams.Harmonization.WaveformPreTrough          = 1.0;   % ms before trough to include
            defaultParams.Harmonization.WaveformPostTrough         = 1.0;   % ms after  trough to include
            defaultParams.Harmonization.WaveformEdgeMode           = "trim"; % "zero" | "edge" | "trim"
                % zero — zero-pad regions outside available input
                % edge — hold the boundary value (nearest-neighbor extrapolation)
                % trim — shrink output window to the shortest available span across all units
            defaultParams.Harmonization.ACGBinSize                 = 0.0005; % s — ACG bin width (0.5 ms)
            defaultParams.Harmonization.ACGLag                     = 0.1;    % s — one-sided ACG lag (100 ms → 401 bins)
            defaultParams.Harmonization.ACGSource                  = "FullACG"; % "FullACG" or "ACG" — which DeePhys ACG property to use
>>>>>>> claude/dreamy-varahamihira
        end

        function acg = computeACG(spike_times, bin_size, lag)
        % COMPUTEACG  Autocorrelogram for a single unit from spike times.
        %
        %   acg = CellTypeClassifier.computeACG(spike_times)
        %   acg = CellTypeClassifier.computeACG(spike_times, bin_size, lag)
        %
        %   spike_times - (N x 1) spike times in seconds
        %   bin_size    - bin width in seconds (default: 0.0005)
        %   lag         - half-width of the ACG window in seconds (default: 0.1)
        %
        %   Returns acg as (N_bins x 1).  Use computeACGs for batches.
            arguments
                spike_times (:,1) double
                bin_size    (1,1) double = 0.0005
                lag         (1,1) double = 0.1
            end
            n_bins = round(lag / bin_size) * 2 + 1;
            acg    = zeros(n_bins, 1);
            if numel(spike_times) > 1
                [acg_u, ~] = CCG(spike_times, ones(numel(spike_times), 1), ...
                    'binSize', bin_size, 'duration', lag * 2);
                if numel(acg_u) == n_bins
                    acg = acg_u;
                end
            end
        end

        function acgs = computeACGs(spike_times_cell, bin_size, lag)
        % COMPUTEACGS  Autocorrelograms for multiple units from spike times.
        %
        %   acgs = CellTypeClassifier.computeACGs(spike_times_cell)
        %   acgs = CellTypeClassifier.computeACGs(spike_times_cell, bin_size, lag)
        %
        %   spike_times_cell - (1 x N) cell array of spike-time vectors (seconds)
        %   bin_size         - bin width in seconds (default: 0.0005)
        %   lag              - half-width of the ACG window in seconds (default: 0.1)
        %
        %   Returns acgs as (N_bins x N), ready for classifyExternalUnits.
            arguments
                spike_times_cell (1,:) cell
                bin_size         (1,1) double = 0.0005
                lag              (1,1) double = 0.1
            end
            N    = numel(spike_times_cell);
            acgs = zeros(round(lag / bin_size) * 2 + 1, N);
            for u = 1:N
                acgs(:, u) = CellTypeClassifier.computeACG( ...
                    spike_times_cell{u}(:), bin_size, lag);
            end
        end

        function stats = evaluateLabels(y_pred, y_true)
        % EVALUATELABELS  Classification performance against ground-truth labels.
        %
        %   stats = CellTypeClassifier.evaluateLabels(y_pred, y_true)
        %
        %   y_pred - (1 x N) predicted labels (1=exc, 2=inh, NaN=unclassified)
        %   y_true - (1 x N) ground-truth labels (1=exc, 2=inh)
        %
        %   NaN entries in either vector are excluded.
        %   Returns a struct with:
        %     accuracy         - overall fraction correct
        %     confusionMatrix  - (2 x 2) rows=true class, cols=predicted class
        %     perClassAccuracy - [acc_exc, acc_inh]
        %     classLabels      - [1, 2]
        %     nValid           - number of non-NaN pairs evaluated
        %     nExcluded        - number of pairs excluded (NaN in either vector)
            arguments
                y_pred (1,:) double
                y_true (1,:) double
            end
            assert(numel(y_pred) == numel(y_true), ...
                'y_pred (%d) and y_true (%d) must have the same length.', ...
                numel(y_pred), numel(y_true));

            valid = ~isnan(y_pred) & ~isnan(y_true);
            yp    = y_pred(valid);
            yt    = y_true(valid);

            classes = [1, 2];
            n_cls   = numel(classes);
            cm      = zeros(n_cls);
            for i = 1:n_cls
                for j = 1:n_cls
                    cm(i, j) = sum(yt == classes(i) & yp == classes(j));
                end
            end

            per_class_acc = zeros(1, n_cls);
            for i = 1:n_cls
                n_true = sum(yt == classes(i));
                if n_true > 0
                    per_class_acc(i) = cm(i, i) / n_true;
                else
                    per_class_acc(i) = NaN;
                end
            end

            stats.accuracy         = sum(diag(cm)) / max(sum(cm(:)), 1);
            stats.confusionMatrix  = cm;
            stats.perClassAccuracy = per_class_acc;
            stats.classLabels      = classes;
            stats.nValid           = sum(valid);
            stats.nExcluded        = sum(~valid);

            fprintf('Accuracy: %.1f%%  |  Exc: %.1f%%  |  Inh: %.1f%%  |  n=%d (excluded: %d NaN)\n', ...
                stats.accuracy * 100, per_class_acc(1) * 100, per_class_acc(2) * 100, ...
                stats.nValid, stats.nExcluded);
        end

    end

    methods (Static, Access = private)

        function [X_tr, X_te, feat_names] = normalizeForClassification(X_tr, X_te, feat_names)
        % NORMALIZEFORCLASSIFICATION  Steps 2–4 of the shared normalization pipeline.
        %
        %   [X_tr, X_te, feat_names] = normalizeForClassification(X_tr, X_te, feat_names)
        %
        %   Applies in sequence:
        %     1. Global z-score: fit on X_tr, apply train statistics to X_te
        %     2. Drop columns with NaN or Inf in either X_tr or X_te
        %     3. Scale by max(abs(X_tr)) column-wise
        %
        %   Per-group z-score (step 1 of the full pipeline) is the caller's
        %   responsibility, as it differs between DeePhys-only and transfer paths.
            [X_tr, mu, sd] = normalize(X_tr, 1, 'zscore');
            X_te = normalize(X_te, 'center', mu, 'scale', sd);

            bad_cols = any(isnan(X_tr) | isinf(X_tr), 1) | ...
                       any(isnan(X_te) | isinf(X_te), 1);
            X_tr(:, bad_cols)    = [];
            X_te(:, bad_cols)    = [];
            feat_names(bad_cols) = [];

            scale = max(abs(X_tr), [], 1);
            scale(scale == 0) = 1;
            X_tr = X_tr ./ scale;
            X_te = X_te ./ scale;
        end

        function mask = buildCultureSubsetMask(fs, culture_indices, culture_keys)
        % BUILDCULTURESUBSETMASK  Logical mask over UnitTable for specified culture indices.
        %   culture_indices are 1-based indices into the list of unique cultures.
        %   culture_keys: string array of metadata fields that define a culture identity
        %   (default ["ChipID","PlatingDate"]).
            arguments
                fs              FeatureStore
                culture_indices (1,:) double
                culture_keys    (1,:) string = ["ChipID", "PlatingDate"]
            end
            meta      = fs.MetadataTable;
            rec_ids   = fs.UnitTable.RecordingID;
            cult_ids  = FeatureStore.getCultureIDsForUnits(rec_ids, meta, culture_keys);
            unique_c  = unique(cult_ids, 'stable');
            mask      = false(1, numel(cult_ids));
            for ci = 1:numel(culture_indices)
                idx = culture_indices(ci);
                if idx >= 1 && idx <= numel(unique_c)
                    mask(cult_ids == unique_c(idx)) = true;
                end
            end
        end

        function ud_array = unitDataFromRecordingGroup(rg)
        % Extract UnitData array from old RecordingGroup (for fromLegacyGroup).
            unit_obj = [rg.Units];
            ud_array = UnitData.empty();
            for u = 1:numel(unit_obj)
                rec_id = string(paramHash(struct('InputPath', ...
                    char(unit_obj(u).MEArecording.Metadata.InputPath))));
                ud_array(u) = UnitData.fromLegacyUnit(unit_obj(u), rec_id);
            end
        end

    end
end
