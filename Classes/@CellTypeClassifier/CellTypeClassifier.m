classdef CellTypeClassifier < handle
% CELLTYPECLASSIFIER  Transductive graph-based cell-type classification pipeline for MEA experiments.
%
% Identifies inhibitory (interneuron) vs excitatory neurons using unsupervised
% unit graph with Leiden / Louvain community detection and graph label propagation.
% Two ground truth strategies are supported:
%
%   Drug response (GroundTruthMethod = "two_window" or "full_curve"):
%     Units that increase firing rate in response to a stimulus serve as
%     inhibitory ground truth. Excitatory counterexamples are inferred from
%     the non-responsive pool via farthest-point sampling.
%
%   Pure culture (GroundTruthMethod = "metadata"):
%     Cell type labels come from a UnitTable column (e.g. CellType field
%     set from culture metadata). Both classes have explicit ground truth.
%     No firing-rate testing is performed.
%
% USAGE:
%   ctc = CellTypeClassifier(featureStore, unitDataArray, params);
%   ctc.identifyResponsiveUnits();    % assign ground truth (FR test or metadata)
%   ctc.generateTrainLabels();        % UMAP embedding + training label assembly
%   ctc.classifyUnits();              % graph label propagation on UMAP graph (transductive)
%   labels = ctc.UnitLabels;          % 1 = excitatory, 2 = inhibitory, NaN = unclassified
%
% MIGRATION from old API:
%   ctc = CellTypeClassifier.fromLegacyGroup(rg, params);

    properties
        FeatureStore            FeatureStore    % Features + metadata for all units
        UnitDataArray           UnitData        % (1 x N) raw unit data (spike times etc.)
        OriginalUnitDataArray   UnitData        % Pre-segment-filter UnitData (all recordings, for FullACG recomputation)
        Parameters              struct          % Merged from returnDefaultParams + user overrides
        ResponsiveUnitIdx       logical         % (1 x N) units with a significant firing rate response
        ResponsiveUnitDirection string          % (1 x N) "none" | "increase" | "decrease" per unit
        ResponsiveStrength      double          % (1 x N) continuous response strength (rho or effect size)
        ResponsivenessDetail    struct          % Per-unit dose-response data for diagnostics
                                               %   .fr_matrix  (n_unique_units x n_doses) firing rates
                                               %   .dose_values (1 x n_doses) GroupingVar values
                                               %   .unit_ids   (1 x n_unique_units) UnitID strings
        CounterexampleUnitIdx   logical         % (1 x N) explicit excitatory ground truth (metadata method only)
        TrainLabels             struct          % .sorted_train_ids, .sorted_y_train, etc.
                                               %   Diagnostic fields: .resp_outlier_mask, .resp_geom_mask,
                                               %   .ce_outlier_mask, .ce_distances_from_inh_centroid,
                                               %   .inh_centroid, .filter_counts
        UnitLabels              double          % (1 x N): 1=excitatory, 2=inhibitory, NaN=unclassified
        UnitConfidence          double          % (1 x N): kNN confidence in [0,1]; 1.0 for training units
        UnitGraphConnectivity   double          % (1 x N): number of direct training-unit neighbors in UMAP graph
        Composition             struct = struct()   % estimateComposition: .Table per culture, .Probability /
                                               %   .Labels per unit (calibrated, culture-prior adjusted)
        UMAP                                   % Trained UMAP model from run_umap
        Graph                   struct = struct()   % Unit graph used for communities and label propagation
                                               %   .W (n_unique x n_unique sparse), .k, .Method ("exact" | "umap")
        Reduction               struct = struct()   % Unsupervised/Train/Test/External embeddings (a 0x0 struct cannot take new fields)
        HarmonizedWaveforms     double          % (N_samples x N_units) processed waveforms
        HarmonizedACGs          double          % (N_bins x N_units) ACGs
        HarmonizedSR            double          % Waveform sampling rate after harmonization
        NormalizationParams     struct          % Normalization pipeline params from generateTrainLabels
                                               %   .mu_global       (1 x F) mean from global z-score
                                               %   .sigma_global    (1 x F) std from global z-score
                                               %   .nan_cols        (1 x F) logical mask of removed NaN columns
                                               %   .scale           (1 x F') max(abs) scaling vector
                                               %   .feat_names_trimmed (1 x F') feature names after NaN removal
        CachedExtraction        struct          % Cached waveform/ACG extraction
        NormalizedFeatures      struct          % Shared preprocessing cache (populated by buildNormalizedFeatures)
                                               %   .X_pergroup   (n_unique x F) after NaN-fill + per-group z-score
                                               %   .feat_names   (1 x F) feature name strings
                                               %   .unique_ud    (1 x n_unique) UnitData, one per unique UnitID
                                               %   .all_to_unique(1 x N_full) index from full array into unique_ud
                                               %   .unique_to_rep(1 x n_unique) representative row in UnitDataArray
                                               %   .subset_mask  (1 x n_unique) training-culture logical mask
                                               %   .n_unique, .n_units_full
    end

    % Kept for backward compatibility (populated by fromLegacyGroup)
    properties (Access = private)
        LegacyRecordingGroup    % RecordingGroup, only set by fromLegacyGroup
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
            ctc.validateParameters();
        end

        function s = saveobj(ctc)
        % SAVEOBJ  Custom serialization — strips non-serializable handles.
        %
        % The UMAP model (ctc.UMAP) is a Java/handle object from run_umap that
        % does not survive save/load across MATLAB sessions. It is stripped here
        % and set to [] on reload. Re-run classifyUnits() after loading if you
        % need a fresh embedding.
        %
        % Computation caches (NormalizedFeatures, CachedExtraction,
        % HarmonizedWaveforms/ACGs) are also stripped — they are regenerated
        % automatically on the next pipeline call.
        %
        % All inference-critical state (TrainLabels, NormalizationParams,
        % UnitLabels, UnitConfidence, FeatureStore, UnitDataArray, Parameters)
        % is preserved.
            s = struct();
            mc = ?CellTypeClassifier;
            for pi = 1:numel(mc.PropertyList)
                pname = mc.PropertyList(pi).Name;
                if mc.PropertyList(pi).GetAccess == "public"
                    s.(pname) = ctc.(pname);
                end
            end
            % Strip the non-serializable UMAP handle and regeneratable caches.
            s.UMAP              = [];
            s.NormalizedFeatures = [];
            s.CachedExtraction  = [];
            s.HarmonizedWaveforms = [];
            s.HarmonizedACGs    = [];
            s.HarmonizedSR      = [];
        end

        function G = graphMatrix(ctc)
        % GRAPHMATRIX  Symmetric unit graph used for communities and label propagation.
        %
        % The exact kNN graph (Parameters.Graph.Method = "exact", deterministic) when it was
        % built, otherwise the fuzzy graph of the UMAP model (randomised neighbour search).
            if isfield(ctc.Graph, 'W') && ~isempty(ctc.Graph.W)
                G = ctc.Graph.W;
            else
                assert(~isempty(ctc.UMAP) && ~isempty(ctc.UMAP.head), ...
                    'CellTypeClassifier:noGraph', ...
                    'No unit graph — run generateTrainLabels() before classification.');
                n = max([ctc.UMAP.head(:); ctc.UMAP.tail(:)]);
                if ~isempty(ctc.Reduction) && isfield(ctc.Reduction, 'Unsupervised') && ~isempty(ctc.Reduction.Unsupervised)
                    n = size(ctc.Reduction.Unsupervised, 1);
                end
                G = sparse(ctc.UMAP.head, ctc.UMAP.tail, ctc.UMAP.weights, n, n);
                G = (G + G') / 2;
            end
        end

        function clearCache(ctc)
        % CLEARCACHE  Invalidate all cached preprocessing and classification state.
        %
        % Call after changing Parameters, ResponsiveUnitIdx, or UnitDataArray
        % post-construction to ensure subsequent pipeline calls recompute from scratch.
        % generateTrainLabels() calls this automatically at its start.
            ctc.NormalizedFeatures  = [];
            ctc.CachedExtraction    = [];
            ctc.HarmonizedWaveforms = [];
            ctc.HarmonizedACGs      = [];
            ctc.TrainLabels         = [];
            ctc.NormalizationParams = [];
            ctc.UnitLabels          = [];
            ctc.UnitConfidence      = [];
            ctc.UnitGraphConnectivity = [];
        end

        function validateParameters(ctc)
        % VALIDATEPARAMETERS  Check for common parameter misconfiguration and warn.
        %
        % Called automatically at the end of the constructor. Can also be called
        % manually after updating ctc.Parameters.
            p      = ctc.Parameters;
            fs     = ctc.FeatureStore;
            n_warn = 0;

            % ── ACG source ───────────────────────────────────────────────────
            if ~isempty(fs) && ~isempty(fs.UnitTable) && ...
                    string(p.Harmonization.ACGSource) == "FullACG"
                all_cols = string(fs.UnitTable.Properties.VariableNames);
                has_acg  = any(startsWith(all_cols, "Parent_ACG") | ...
                               startsWith(all_cols, "FullACG"));
                if ~has_acg
                    warning('CellTypeClassifier:noFullACGColumns', ...
                        ['Harmonization.ACGSource="FullACG" but no Parent_ACG* or ' ...
                         'FullACG* columns found in FeatureStore.UnitTable. ' ...
                         'Will fall back to segment-level ACG computation.']);
                    n_warn = n_warn + 1;
                end
            end

            % ── Recordings per culture for full_curve ────────────────────────
            if string(p.Bootstrap.GroundTruthMethod) == "full_curve" && ...
                    ~isempty(fs) && ~isempty(fs.MetadataTable)
                meta     = fs.MetadataTable;
                cult_ids = FeatureStore.getCultureIDsForUnits( ...
                    fs.UnitTable.RecordingID, meta, p.CultureKeys);
                unique_c = unique(cult_ids, 'stable');
                n_recs_per_cult = arrayfun(@(c) ...
                    numel(unique(fs.UnitTable.RecordingID(cult_ids == c))), unique_c);
                n_few = sum(n_recs_per_cult < p.Bootstrap.MinRecordings);
                if n_few > 0
                    warning('CellTypeClassifier:fewRecordingsForFullCurve', ...
                        ['%d/%d cultures have fewer than MinRecordings=%d recordings. ' ...
                         'They will fall back to the two_window method.'], ...
                        n_few, numel(unique_c), p.Bootstrap.MinRecordings);
                    n_warn = n_warn + 1;
                end
            end

            % ── TemplateDir ─────────────────────────────────────────────────
            td = p.UMAP.TemplateDir;
            if ~isempty(td) && ~isfolder(td)
                warning('CellTypeClassifier:templateDirMissing', ...
                    ['UMAP.TemplateDir "%s" does not exist. ' ...
                     'Create the directory or set Parameters.UMAP.TemplateDir to a valid path.'], td);
                n_warn = n_warn + 1;
            end

            % ── NNeighbors vs dataset size ───────────────────────────────────
            if ~isempty(fs) && ~isempty(fs.UnitTable)
                n_unique = numel(unique(string(fs.UnitTable.UnitID)));
                if p.UMAP.NNeighbors >= n_unique
                    warning('CellTypeClassifier:nNeighborsTooLarge', ...
                        'UMAP.NNeighbors (%d) >= number of unique units (%d). Will be clamped.', ...
                        p.UMAP.NNeighbors, n_unique);
                    n_warn = n_warn + 1;
                end
            end

            if n_warn == 0
                % No issues found — silent by design.
            end
        end

        function [wf, acg, sr] = getOrExtract(ctc, unit_data)
        % GETOREXTRACT  Three-level cache: in-memory -> disk -> compute.
        %
        % ACG source strategy (ACGSource parameter):
        %   "FullACG"      — prefer Parent_ACG* columns from FeatureStore (pre-computed
        %                    from full parent spike train). Falls back to the
        %                    concatenated segments, then to the segment ACG.
        %   "Concatenated" — spike trains of all segments of the unit joined.
        %   "ACG"          — the unit's own segment only (UnitData.SpikeTimes). For
        %                    dose-response data this is the drug-free (baseline) ACG;
        %                    the other two include the drug periods.
            ph    = ctc.Parameters.Harmonization;
            N     = numel(unit_data);
            cache = ctc.CachedExtraction;

            % -- Level 1: in-memory cache -----------------------------------------
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

            % -- Level 2: disk cache ----------------------------------------------
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
                % -- Level 3: compute ---------------------------------------------
                % Waveform: always from UnitData.
                % Units from different recordings may have different waveform lengths
                % (different sampling rates or pre/post trough cutouts). Collect
                % individually and pad shorter waveforms with NaN so buildFeatureMatrix
                % can resample everything to the common target grid.
                wf_list   = {unit_data.ReferenceWaveform};
                wf_lens   = cellfun(@numel, wf_list);
                max_len   = max(wf_lens);
                if max_len == 0
                    wf = zeros(0, N);
                    sr = ph.WaveformTargetSamplingRate;
                else
                    wf = nan(max_len, N);
                    for i_wf = 1:N
                        L = wf_lens(i_wf);
                        if L > 0
                            wf(1:L, i_wf) = double(wf_list{i_wf}(:));
                        end
                    end
                    sr = [unit_data.SamplingRate];   % per unit: datasets may mix rates
                end

                % ACG: prefer Parent_ACG* from FeatureStore when ACGSource == "FullACG"
                acg = ctc.resolveACG(unit_data, ph);

                % Write to disk cache
                try
                    save(cache_file, 'wf', 'acg', 'sr', '-v7.3');
                catch
                end
            end

            % -- Update in-memory cache -------------------------------------------
            ctc.CachedExtraction = struct('wf', wf, 'acg', acg, 'sr', sr, ...
                'N', N, 'ACGBinSize', ph.ACGBinSize, 'ACGLag', ph.ACGLag, ...
                'ACGSource', ph.ACGSource);
        end

        function acg = resolveACG(ctc, unit_data, ph)
        % RESOLVEACG  Return ACG matrix, preferring full-recording ACG columns from FeatureStore.
        %
        % Accepts either "FullACG*" (legacy MEArecording.mat convention) or
        % "Parent_ACG*" column names — whichever is present in the UnitTable.
            if ph.ACGSource == "FullACG" && ~isempty(ctc.FeatureStore)
                ut       = ctc.FeatureStore.UnitTable;
                all_cols = string(ut.Properties.VariableNames);
                % Accept both legacy "FullACG*" and newer "Parent_ACG*" naming
                parent_acg_cols = all_cols(startsWith(all_cols, "Parent_ACG") | ...
                                           startsWith(all_cols, "FullACG"));
                if ~isempty(parent_acg_cols)
                    % Map unit_data order to FeatureStore row order by (UnitID, RecordingID)
                    % compound key. UnitID alone is ambiguous in dose-response data where
                    % the same unit appears in multiple recordings.
                    ud_keys   = string({unit_data.UnitID})' + "|" + string({unit_data.RecordingID})';
                    fs_keys   = string(ut.UnitID) + "|" + string(ut.RecordingID);
                    [~, loc]  = ismember(ud_keys, fs_keys);
                    valid     = loc > 0;
                    n_bins    = numel(parent_acg_cols);
                    acg       = nan(n_bins, numel(unit_data));
                    if any(valid)
                        acg(:, valid) = ut{loc(valid), parent_acg_cols}';
                    end
                    if any(~valid)
                        warning('CellTypeClassifier:acgMissing', ...
                            '%d/%d units have no Parent_ACG row in the FeatureStore; their ACG is NaN.', ...
                            sum(~valid), numel(valid));
                    end
                    % Stored ACGs use the processor's window; adapt or recompute when it differs.
                    n_target = round(2 * ph.ACGLag / ph.ACGBinSize) + 1;
                    if n_bins == n_target
                        return
                    end
                    [from_bs, from_lag] = CellTypeClassifier.storedACGParams(unit_data, n_bins);
                    adapted = CellTypeClassifier.adaptACG(acg, from_bs, from_lag, ph.ACGBinSize, ph.ACGLag);
                    if ~isempty(adapted)
                        fprintf('CellTypeClassifier: Parent_ACG adapted from %d bins (%.4g s, +-%.4g s) to %d bins.\n', ...
                            n_bins, from_bs, from_lag, n_target);
                        acg = adapted;
                        return
                    end
                    fprintf(['CellTypeClassifier: Parent_ACG (%d bins) cannot be adapted to Harmonization ' ...
                        '(BinSize=%.4g, Lag=%.4g) — recomputing from spike trains.\n'], n_bins, ph.ACGBinSize, ph.ACGLag);
                end
                % FeatureStore has no FullACG columns — try reading FullACG
                % directly from the UnitData structs (legacy MEArecording path).
                has_fullacg = ~isempty(unit_data) && ...
                    ((isstruct(unit_data) && isfield(unit_data, 'FullACG')) || ...
                     (isobject(unit_data) && isprop(unit_data(1), 'FullACG'))) && ...
                    ~isempty(unit_data(1).FullACG);
                if has_fullacg
                    % Check whether stored FullACG parameters match requested Harmonization.
                    % Compare bin size and lag explicitly (not just bin count).
                    stored_bs  = unit_data(1).FullACGBinSize;
                    stored_lag = unit_data(1).FullACGLag;
                    params_known = ~isempty(stored_bs) && ~isempty(stored_lag);
                    if params_known
                        params_match = abs(stored_bs - ph.ACGBinSize) < 1e-9 ...
                                    && abs(stored_lag - ph.ACGLag) < 1e-9;
                    else
                        % Parameters not stored — fall back to bin count check
                        n_bins_target = round(2 * ph.ACGLag / ph.ACGBinSize) + 1;
                        params_match  = numel(unit_data(1).FullACG) == n_bins_target;
                    end

                    if params_match
                        n_bins = round(2 * ph.ACGLag / ph.ACGBinSize) + 1;
                        acg = zeros(n_bins, numel(unit_data));
                        n_recomputed = 0;
                        for u = 1:numel(unit_data)
                            fa = unit_data(u).FullACG;
                            if ~isempty(fa) && numel(fa) == n_bins
                                acg(:, u) = fa(:);
                            else
                                % Stored FullACG has wrong size — recompute from parent or segment
                                if ~isempty(ctc.OriginalUnitDataArray) && numel(ctc.OriginalUnitDataArray) > numel(unit_data)
                                    uid = string(unit_data(u).UnitID);
                                    seg_idx = find(string({ctc.OriginalUnitDataArray.UnitID}) == uid);
                                    combined = [];
                                    offset = 0;
                                    gap = 2 * ph.ACGLag;
                                    for si = 1:numel(seg_idx)
                                        st = ctc.OriginalUnitDataArray(seg_idx(si)).SpikeTimes(:);
                                        if ~isempty(st)
                                            combined = [combined; st + offset]; %#ok<AGROW>
                                        end
                                        offset = offset + ctc.OriginalUnitDataArray(seg_idx(si)).RecordingDuration + gap;
                                    end
                                    acg(:, u) = CellTypeClassifier.computeACG(combined, ph.ACGBinSize, ph.ACGLag);
                                else
                                    acg(:, u) = CellTypeClassifier.computeACG( ...
                                        unit_data(u).SpikeTimes, ph.ACGBinSize, ph.ACGLag);
                                end
                                n_recomputed = n_recomputed + 1;
                            end
                        end
                        if n_recomputed > 0
                            fprintf('CellTypeClassifier: FullACG from UnitData (%d bins), %d/%d recomputed (size mismatch).\n', ...
                                n_bins, n_recomputed, numel(unit_data));
                        else
                            fprintf('CellTypeClassifier: FullACG from UnitData (%d bins, BinSize=%.4f, Lag=%.4f).\n', ...
                                n_bins, ph.ACGBinSize, ph.ACGLag);
                        end
                        return
                    else
                        if params_known
                            fprintf(['CellTypeClassifier: FullACG params (BinSize=%.4f, Lag=%.4f) ' ...
                                'differ from Harmonization (BinSize=%.4f, Lag=%.4f) — recomputing.\n'], ...
                                stored_bs, stored_lag, ph.ACGBinSize, ph.ACGLag);
                        else
                            fprintf(['CellTypeClassifier: FullACG bin count (%d) does not match ' ...
                                'Harmonization (%d) and no stored params — recomputing.\n'], ...
                                numel(unit_data(1).FullACG), round(2*ph.ACGLag/ph.ACGBinSize)+1);
                        end
                    end
                else
                    fprintf('CellTypeClassifier: No FullACG in FeatureStore or UnitData — using segment ACG.\n');
                end
            end
            % "FullACG" without stored ACGs, and "Concatenated": join the spike trains of all
            % segments (OriginalUnitDataArray, kept by discardSegmentUnits) to approximate the
            % parent recording. More spikes, but in dose-response data this includes the drug
            % periods, so the ACG of a responsive unit carries its drug response.
            % "ACG": the unit's own (baseline) segment only, i.e. drug-free features.
            use_concat = ismember(string(ph.ACGSource), ["FullACG", "Concatenated"]);
            if use_concat && ~isempty(ctc.OriginalUnitDataArray) && numel(ctc.OriginalUnitDataArray) > numel(unit_data)
                acg = ctc.computeParentACGs(unit_data, ph);
                fprintf('CellTypeClassifier: ACG computed from parent spike trains (%d units, BinSize=%.4f, Lag=%.4f).\n', ...
                    numel(unit_data), ph.ACGBinSize, ph.ACGLag);
            else
                % Fall back to segment-level computation from SpikeTimes
                [~, acg, ~] = ctc.extractFromUnitData(unit_data);
            end
        end

    end

    % =====================================================================
    % Static constructors and defaults
    % =====================================================================
    methods (Static)

        function ctc = loadobj(s)
        % LOADOBJ  Reconstruct CellTypeClassifier from a saved struct.
        %
        % Called automatically by MATLAB load(). Restores all serialized
        % properties from the struct written by saveobj(), then warns if
        % UnitDataArray is empty (required by classifyExternalUnits and
        % classifyUnitsWithExternalTrain).
            ctc = CellTypeClassifier();
            if isstruct(s)
                fn = fieldnames(s);
                for fi = 1:numel(fn)
                    try
                        ctc.(fn{fi}) = s.(fn{fi});
                    catch
                        % Skip properties that no longer exist in the class
                    end
                end
            end
            if isempty(ctc.UnitDataArray)
                warning('CellTypeClassifier:loadobj', ...
                    ['UnitDataArray is empty after load. ' ...
                     'classifyExternalUnits() and classifyUnitsWithExternalTrain() ' ...
                     'require UnitDataArray — restore it via ctc.UnitDataArray = ud ' ...
                     'before calling those methods.']);
            end
        end

        function ctc = fromProcessors(proc_array, parameters)
        % FROMPROCESSORS  Build CellTypeClassifier from RecordingProcessor array.
        %   ctc = CellTypeClassifier.fromProcessors(procs, params)
            arguments
                proc_array  RecordingProcessor
                parameters  struct = struct()
            end
            fs  = FeatureStore.fromProcessors(proc_array);
            % Match FeatureStore.fromProcessors inclusion logic exactly: skip processors
            % with empty SpikeData (which FeatureStore skips entirely) and processors
            % with no features/units (which contribute zero unit rows).
            has_features = arrayfun(@(p) ...
                ~isempty(p.SpikeData) && ~isempty(p.SpikeData.RecordingID) && ...
                ~isempty(p.UnitFeatureTable) && ~isempty(p.Units), ...
                proc_array);
            ud  = [proc_array(has_features).Units];
            ctc = CellTypeClassifier(fs, ud, parameters);
        end

        function ctc = fromDatabase(parameters, opts)
        % FROMDATABASE  Query RecordingDatabase and build CellTypeClassifier.
        %   ctc = CellTypeClassifier.fromDatabase(params, Mutation="WT", DIV=[14 21])
        %
        %   Finds processor files at <input_path>/RecordingProcessor.mat or
        %   <input_path>/MEArecording.mat (legacy). RecordingProcessor.load()
        %   auto-migrates legacy format.
            arguments
                parameters  struct = struct()
                opts.Mutation string = ""
                opts.ChipID string = ""
                opts.PlateID string = ""
                opts.WellID string = ""
                opts.DIV double = []
            end
            db = RecordingDatabase.instance();
            query_args = {};
            if opts.Mutation ~= ""; query_args = [query_args, {'Mutation', opts.Mutation}]; end
            if opts.ChipID ~= "";   query_args = [query_args, {'ChipID', opts.ChipID}]; end
            if opts.PlateID ~= "";  query_args = [query_args, {'PlateID', opts.PlateID}]; end
            if opts.WellID ~= "";   query_args = [query_args, {'WellID', opts.WellID}]; end
            if ~isempty(opts.DIV);  query_args = [query_args, {'DIV', opts.DIV}]; end
            T = db.queryRecordings(query_args{:});
            assert(~isempty(T) && height(T) > 0, ...
                'CellTypeClassifier:noRecordings', 'No recordings match the query.');

            input_paths = string(T.input_path);
            proc_paths = strings(numel(input_paths), 1);
            for i = 1:numel(input_paths)
                ip = input_paths(i);
                rp = fullfile(ip, "RecordingProcessor.mat");
                mp = fullfile(ip, "MEArecording.mat");
                if isfile(rp)
                    proc_paths(i) = rp;
                elseif isfile(mp)
                    proc_paths(i) = mp;
                elseif isfile(ip) && endsWith(ip, ".mat")
                    proc_paths(i) = ip;
                end
            end
            found = proc_paths ~= "";
            assert(any(found), 'CellTypeClassifier:noFiles', ...
                'No processor/MEArecording files found for queried recordings.');
            if ~all(found)
                fprintf('fromDatabase: %d/%d files found, skipping %d missing\n', ...
                    sum(found), numel(found), sum(~found));
            end
            procs = RecordingProcessor.loadMany(proc_paths(found));
            ctc = CellTypeClassifier.fromProcessors(procs, parameters);
        end

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
        % RETURNDEFAULTPARAMS  Default parameters. Most of them only serve optional methods.
        %
        % The default pipeline ("full_curve" responders, Graph.Method "exact", "enrichment"
        % seeds, graph label propagation) reads:
        %   Experiment description (must match the data; nothing to tune):
        %     CultureKeys; UMAP.GroupingVar / GroupingValues (dose variable, baseline value);
        %     UMAP.NormalizationVar (per-recording z-score for cultures of natural composition,
        %     "" for a set enriched in one class); Bootstrap.Direction;
        %     TrainLabels.ResponsiveClassLabel; the candidate cultures are the filter passed
        %     to identifyResponsiveUnits.
        %   Feature definition (fixed by what the data contain): Harmonization.* (waveform
        %     window and rate, ACG bins and lag, ACGSource); UMAP.ACGWeight / WaveformWeight.
        %   Chosen from the data: Harmonization.ACGSmoothing ("auto").
        %   Fixed settings, not tuned per dataset:
        %     responders  Bootstrap.MonotonicityThreshold, MinFoldChange, ConfirmationAlpha,
        %                 ConfirmationNIter, NetworkCorrection, MinRecordings
        %     seeds       UMAP.NNeighbors (graph k), Community.LouvainResolution,
        %                 LouvainRestarts, SeedFDR
        % Everything else belongs to alternative or legacy paths: the two_window /
        % control_calibrated / metadata responder settings, the UMAP embedding (NDims, MinDist,
        % Spread: plots only with the exact graph), OutlierDetection.* and the remaining
        % Community.* settings of the "community" / "iforest" seed methods, the "knn"
        % classifier, Ensemble, BayesianOptimization and the supervised-UMAP settings.
        % qualityReport() replaces parameter tuning by acceptance checks.

            % ── Bootstrap: responsive unit identification ─────────────────────
            % GroundTruthMethod: how to identify inhibitory candidates.
            %   "two_window" — compare pre vs post spike rates within a single recording
            %                  (bootstrap permutation test). Uses Pre/PostCutout, BinSize,
            %                  NIter, Alpha.
            %   "full_curve" — per-unit Spearman rank correlation between firing rate and
            %                  GroupingVar (e.g. Concentration) across all recordings in the
            %                  culture. Same UnitIDs must appear across recordings. Uses
            %                  FullCurveAlpha. Pre/PostCutout, BinSize, NIter, Alpha not used.
            %                  Falls back to "two_window" if fewer than MinRecordings recordings.
            %   "metadata"    — per-unit cell type labels come from a column in FeatureStore.UnitTable.
            %                  LabelField specifies the column; ResponsiveClassValue specifies which
            %                  value maps to ResponsiveClassLabel. Units with other non-empty values
            %                  become explicit counterexample ground truth. No FR tests are run.
            defaultParams.Bootstrap.GroundTruthMethod     = "full_curve";
            % metadata: read cell type labels directly from UnitTable.
            %   LabelField: column name in FeatureStore.UnitTable (e.g. "CellType")
            %   ResponsiveClassValue: value in that column that maps to ResponsiveClassLabel
            %     (e.g. "inhibitory"). Units with other non-empty values become counterexamples.
            defaultParams.Bootstrap.LabelField            = "";   % e.g. "CellType"
            defaultParams.Bootstrap.ResponsiveClassValue  = "";   % e.g. "inhibitory"
            % control_calibrated: per unit, the firing-rate change from baseline to the
            %   HighDoses highest doses (log2 fold change minus the culture median, rank-sum
            %   test on BinSize-s counts). Thresholds are the ControlQuantile tails of the same
            %   quantity in control cultures (ControlFilter, e.g. {'Construct', 0}: no DREADD),
            %   so the call rate in controls is the false-positive rate and the estimated FDR is
            %   reported. Replaces Alpha/NIter/NetworkCorrection/MinFoldChange tuning; needs
            %   control cultures recorded with the same protocol.
            defaultParams.Bootstrap.ControlFilter         = {};   % e.g. {'Construct', 0}
            defaultParams.Bootstrap.ControlQuantile       = 0.01;
            defaultParams.Bootstrap.CalibratedAlpha       = 0.01;
            defaultParams.Bootstrap.HighDoses             = 2;
            % TargetFDR: when set, the threshold is the least strict one whose control-estimated
            %   FDR stays <= TargetFDR (more responders at a stated purity); NaN = fixed
            %   ControlQuantile tails. Strict tails keep only the strongest responders, which
            %   can be a biased subset of the target cells.
            defaultParams.Bootstrap.TargetFDR             = NaN;
            % two_window: compares mean FR between two recordings within a culture,
            %   selected by their GroupingVar value.
            %   PreGroupValue/PostGroupValue identify which recordings to compare.
            %   [] = auto-select (lowest GroupingVar for pre, highest for post).
            %   PreCutout/PostCutout optionally restrict to a sub-window (seconds)
            %   within each recording; [] = full recording duration.
            defaultParams.Bootstrap.PreGroupValue         = [];   % [] = lowest GroupingVar recording
            defaultParams.Bootstrap.PostGroupValue        = [];   % [] = highest GroupingVar recording
            defaultParams.Bootstrap.PreCutout             = [];   % [] = full recording; or [start, end] in seconds
            defaultParams.Bootstrap.PostCutout            = [];   % [] = full recording; or [start, end] in seconds
            defaultParams.Bootstrap.MinRecordings         = 4;
            defaultParams.Bootstrap.BinSize               = 20;           % two_window only (seconds)
            defaultParams.Bootstrap.NIter                 = 1000;         % two_window only
            defaultParams.Bootstrap.Alpha                 = 1e-3;         % two_window only; 1e-10 was unreliable for 1000 bootstrap samples
            defaultParams.Bootstrap.MonotonicityThreshold = 0.75;        % full_curve stage 1: min |Spearman rho|
            defaultParams.Bootstrap.MinFoldChange         = 1.5;         % full_curve stage 2: min fold change (baseline to max dose)
            defaultParams.Bootstrap.ConfirmationAlpha     = 1e-3;        % full_curve stage 3: bootstrap p-value threshold
            defaultParams.Bootstrap.ConfirmationNIter     = 1000;        % full_curve stage 3: bootstrap iterations
            defaultParams.Bootstrap.Direction             = "increase";   % "increase" | "decrease" | "both"
            % ResponseRules: several perturbations of the responsive class, one row
            %   {field, value, direction} each, e.g. {'Construct', 1, "increase"; 'Construct', 2, "decrease"}
            %   (activating and silencing DREADD in interneurons). Responders of all rules are
            %   pooled; validateResponseRules trains on each rule alone and validates on the
            %   others - pool only if both directions hold; otherwise keep one rule as the
            %   independent check.
            %   Empty = one condition (identifyResponsiveUnits' metadata_filter + Direction).
            defaultParams.Bootstrap.ResponseRules         = {};
            defaultParams.Bootstrap.MinResponsiveStrength = 0;  % 0 = no filtering

            % UseFDR: replace fixed Alpha with Benjamini-Hochberg FDR control.
            %   Scales correctly with dataset size — large datasets get proper
            %   multiple-comparison correction while small datasets retain power.
            defaultParams.Bootstrap.UseFDR               = false;
            defaultParams.Bootstrap.FDRLevel             = 0.05;

            % UnitWisePermutation: use independent randperm per unit (slower but correct
            %   for correlated spike trains; avoids inflated FPR from shared null).
            defaultParams.Bootstrap.UnitWisePermutation  = true;

            % NetworkCorrection: subtract population-mean FR change before per-unit testing.
            %   Removes the network-wide disinhibition confound (e.g. GABAergic blockade
            %   drives all units to fire more; only units firing above-average are flagged).
            defaultParams.Bootstrap.NetworkCorrection     = true;

            % ── UMAP parameters ───────────────────────────────────────────────
            defaultParams.UMAP.NDims                = 5;
            defaultParams.UMAP.NNeighbors           = 50;
            defaultParams.UMAP.MinDist              = 0.1;
            defaultParams.UMAP.Spread               = 1;
            defaultParams.UMAP.SupervisedNDims      = 2;
            defaultParams.UMAP.SupervisedNNeighbors = 100;
            defaultParams.UMAP.UnitFeatures         = ["FullACG","ReferenceWaveform"];
            defaultParams.UMAP.NormalizationVar     = "ChipID";
            % {field, value}: z-score each NormalizationVar group with the statistics of its
            % reference units only (bridge wells of equal composition on every plate / batch)
            defaultParams.UMAP.NormalizationReference = {};
            defaultParams.UMAP.GroupingVar          = "Concentration";
            defaultParams.UMAP.GroupingValues       = 0;
            defaultParams.UMAP.TrainingCultureIdx   = [];
            defaultParams.UMAP.TargetWeight         = 0.3;  % 0=data topology, 1=label topology
            defaultParams.UMAP.TemplateDir          = '';   % must be set by user; UMAP template saved here

            % ACGWeight / WaveformWeight: dimension-normalised feature group scaling.
            %   Each group's total L2 contribution is proportional to its weight
            %   regardless of group size. sqrt(w/n) is applied column-wise so that
            %   sum-of-squared contributions = w for each group.
            %   1.0 = equal contribution per group (default); 2.0 = double contribution.
            %   Applied before the unsupervised UMAP (runUnsupervisedUMAP and
            %   optimizeUnsupervisedUMAP alike) and before the kNN classifier.
            defaultParams.UMAP.ACGWeight            = 1.0;
            defaultParams.UMAP.WaveformWeight       = 1.0;
            defaultParams.UMAP.ConfidenceK          = 15;  % k for kNN confidence scoring

            % ── Unit graph ───────────────────────────────────────────────────
            % Graph.Method: graph used for communities and label propagation.
            %   "exact" (default) — UMAP's fuzzy kNN graph built with an exact neighbour
            %                       search (buildKNNGraph): identical across runs and seeds.
            %   "umap"            — the graph of the UMAP model (approximate, randomised
            %                       neighbour search; results varied between runs).
            %   "wnn"             — experimental weighted-nearest-neighbour graph (buildWNNGraph;
            %                       Seurat WNN / PhysMAP): ACG and waveform get per-unit
            %                       weights from how well each modality predicts the unit's
            %                       neighbours; independent of UMAP.ACGWeight / WaveformWeight.
            %   With "exact", the UMAP embedding (NDims, MinDist, Spread) only serves plots and
            %   the embedding-based outlier methods ("community", "iforest"); it does not affect
            %   labels from the default "enrichment" method.
            defaultParams.Graph.Method              = "exact";
            % ComputeEmbedding: run UMAP for the 2-D/ND embedding (diagnostic plots). Required
            % by the embedding-based outlier methods; optional for "enrichment".
            defaultParams.UMAP.ComputeEmbedding     = true;

            % AutoNNeighbors: UMAP n_neighbors set to max(MinNNeighbors, sqrt(N)).
            %   Structural scaling: prevents n_neighbors from exceeding dataset size
            %   (small datasets) or being too local (large datasets). Standard heuristic.
            defaultParams.UMAP.AutoNNeighbors       = false;
            defaultParams.UMAP.MinNNeighbors        = 15;

            % AutoConfidenceK: kNN confidence k set to max(5, sqrt(N_train)).
            %   Structural scaling: prevents k from approaching N_train (making all
            %   confidences identical) while scaling with training set density.
            defaultParams.UMAP.AutoConfidenceK      = false;

            % FeatureSelection: remove low-variance and highly correlated features before UMAP.
            %   Data-quality step: reduces ~640-feature input to ~200-400 informative features,
            %   improving UMAP stability across seeds and reducing computation time.
            %   MinVariancePercentile: drop features in the bottom N% by variance.
            %   MaxCorrelation: drop one of each pair with |r| above this threshold.
            defaultParams.UMAP.FeatureSelection      = false;
            defaultParams.UMAP.MinVariancePercentile = 10;
            defaultParams.UMAP.MaxCorrelation        = 0.99;

            % ── Training label generation ─────────────────────────────────────
            defaultParams.TrainLabels.ResponsiveClassLabel = 2;  % 2=responsive->inhibitory, 1=responsive->excitatory

            % ── Outlier detection ─────────────────────────────────────────────
            % Method: controls how responsive-unit outliers are detected and how
            %   counterexamples are selected when no explicit CE labels are provided.
            %   "enrichment" (default) — communities (Community.Algorithm) on the unit graph;
            %                  significantly enriched in responders (hypergeometric test,
            %                  BH-FDR Community.SeedFDR) give the inhibitory seeds, significantly
            %                  depleted ones the excitatory seeds. No thresholds to tune and no
            %                  fallback: errors if no community is enriched (selectEnrichmentSeeds).
            %   "community" — Louvain community filter + graph purity check + k-medoids CE.
            %                  Falls back to "iforest" if communities are not well-separated.
            %   "iforest"   — Isolation forest on feature-space (PCA or UMAP domain) +
            %                  geometric consistency filter + distance-based CE selection.
            %   "none"      — No outlier detection; all responsive units used as training
            %                  labels; distance-based CE selection without iforest filtering.
            defaultParams.OutlierDetection.Method                            = "enrichment";
            defaultParams.OutlierDetection.CounterexampleRatio               = 1;    % 1:1 balanced default; use actual ratio for metadata
            defaultParams.OutlierDetection.CounterexampleDistancePercentile  = 50;   % percentile of inh self-distances used as CE distance threshold
            defaultParams.OutlierDetection.Domain                = "umap";
            defaultParams.OutlierDetection.ContaminationFraction = "auto";
            defaultParams.OutlierDetection.AutoGMMSeparation     = 2;
            defaultParams.OutlierDetection.NPCAComponents        = 15;

            % ── Community detection ───────────────────────────────────────────
            % Used in generateTrainLabels for community-based CE selection.
            %   LouvainResolution: Louvain gamma parameter (1 = standard modularity;
            %     higher = finer communities). Set by optimizeUnsupervisedUMAP or manually.
            %   InhibitoryCommunityRelThresh: community is considered inhibitory if its
            %     responsive fraction >= this × max responsive fraction across communities.
            %   FallbackCriterion: when to abandon the community path for iforest CE selection.
            %     "permutation" (default): fall back unless responsive units concentrate in
            %     the selected communities more than randomly drawn units do (p <=
            %     FallbackAlpha over FallbackPermutations redraws; the partition is fixed).
            %     "absolute": fall back if fewer than CommunityFallbackThreshold of the
            %     responsive units land in inhibitory communities. With many false-positive
            %     candidates that fraction is itself near 0.5, so the path can flip with the seed.
            %   MinCEPerCommunity: minimum k-medoids selected from each CE community.
            %   PuritySigmaThreshold: robust z-score cutoff for graph purity outlier removal.
            %   LouvainRestarts: number of Louvain restarts; best-Q run is used. The UMAP
            %     graph itself does not depend on the seed, so with enough restarts the
            %     best-Q partition (and the labels) become seed-independent.
            %   EnrichmentFactor: community must have >= EnrichmentFactor × p_resp responsive
            %     fraction to qualify as inhibitory (chance-anchored absolute floor).
            defaultParams.Community.LouvainResolution             = 1.0;
            defaultParams.Community.InhibitoryCommunityRelThresh  = 0.3;
            defaultParams.Community.FallbackCriterion             = "permutation";
            defaultParams.Community.FallbackAlpha                 = 0.01;
            defaultParams.Community.FallbackPermutations          = 200;
            defaultParams.Community.CommunityFallbackThreshold    = 0.5;   % "absolute" criterion only
            defaultParams.Community.MinCEPerCommunity             = 3;
            defaultParams.Community.PuritySigmaThreshold          = 2.5;
            defaultParams.Community.LouvainRestarts               = 30;
            defaultParams.Community.EnrichmentFactor              = 1.5;
            % SeedFDR: false-discovery rate for responder enrichment / depletion of communities
            %   ("enrichment" method; replaces EnrichmentFactor, InhibitoryCommunityRelThresh,
            %   PuritySigmaThreshold and the fallback settings above).
            defaultParams.Community.SeedFDR                       = 0.05;
            % Algorithm: community detection for the enrichment seeds. "leiden" (default;
            %   Traag et al. 2019, leidenCommunities: connected communities, works on the
            %   sparse graph) or "louvain" (BCT community_louvain, needs the full N x N
            %   matrix: ~20 GB at 50,000 units). Partition quality is similar; Leiden is
            %   faster and its restarts agree more.
            %   LouvainResolution / LouvainRestarts apply to both.
            defaultParams.Community.Algorithm                     = "leiden";

            % ── Classification ────────────────────────────────────────────────
            defaultParams.Classification.Method                 = "graph";   % "knn" (feature-space kNN) or "graph" (UMAP graph label propagation)
            defaultParams.Classification.GraphMaxIter            = 100;    % graph method: max label propagation iterations
            defaultParams.Classification.GraphConvergenceTol     = 1e-4;   % graph method: convergence tolerance
            defaultParams.Classification.UseConfidenceThreshold  = false;  % mark low-confidence predictions as NaN
            defaultParams.Classification.ConfidenceThreshold     = 0.8;    % threshold when UseConfidenceThreshold=true

            % ── Ensemble classification ───────────────────────────────────────
            % Enabled: run the pipeline N times with different RNG seeds and take
            %   majority-vote labels. Confidence = fraction of seeds agreeing.
            %   MinAgreement: minimum vote fraction for label assignment.
            %   Units below this threshold become NaN (unclassified).
            % Off by default: with Graph.Method "exact" and "enrichment" seeds (a consensus
            % over Community.LouvainRestarts) the labels hardly depend on the seed, so the
            % ensemble only multiplies the runtime and replaces P(inhibitory) by a vote
            % fraction. Useful with Graph.Method "umap" or the embedding-based seed methods.
            defaultParams.Ensemble.Enabled      = false;
            defaultParams.Ensemble.Seeds        = [42, 1042, 2042, 3042, 4042];
            defaultParams.Ensemble.MinAgreement = 0.6;

            % ── Bayesian optimization ────────────────────────────────────────
            % Single-phase optimization (Phase 1):
            %   optimizeUnsupervisedUMAP: optimises unsupervised UMAP + Louvain parameters
            %     for community_coherence — fraction of responsive units in inhibitory
            %     Louvain communities. Variables: UnsupOptimizeVars.
            %   AutoLouvainRange: when true, the LouvainResolution search range is
            %     adapted based on N_all/n_responsive to target communities of
            %     appropriate granularity. User range (LouvainResolutionRange) sets
            %     the hard bounds; AutoLouvainRange adjusts the upper bound.

            defaultParams.BayesianOptimization.MinDistRange              = [0.01, 1.0];
            defaultParams.BayesianOptimization.SpreadRange               = [0.5, 5.0];
            defaultParams.BayesianOptimization.ACGWeightRange            = [0.5, 3.0];
            defaultParams.BayesianOptimization.WaveformWeightRange       = [0.5, 3.0];
            defaultParams.BayesianOptimization.UnsupOptimizeVars      = [ ...
                "NNeighbors", "MinDist", "Spread", "ACGWeight", "WaveformWeight", "LouvainResolution"];
            defaultParams.BayesianOptimization.NNeighborsRange        = [10, 100];
            defaultParams.BayesianOptimization.LouvainResolutionRange = [0.5, 3.0];
            defaultParams.BayesianOptimization.LouvainRestarts        = 5;   % per evaluation (final run: Community.LouvainRestarts)
            defaultParams.BayesianOptimization.AutoLouvainRange       = true;

            % ── Diagnostics ───────────────────────────────────────────────────
            % Enable: when true, each pipeline method generates a diagnostic figure
            %   at the end of its execution. Each diagnostic method is also callable
            %   independently: ctc.diagnosticResponsiveUnits(), etc.
            % SaveDir: when non-empty, figures are saved as PNG to this directory.
            % ShowFigures: set false for headless/batch runs (saves without display).
            defaultParams.Diagnostics.Enable      = false;
            defaultParams.Diagnostics.SaveDir     = '';
            defaultParams.Diagnostics.ShowFigures = true;

            % ── General ───────────────────────────────────────────────────────
            defaultParams.RNGSeed     = 42;
            defaultParams.CultureKeys = ["ChipID", "PlatingDate"];

            % ── Harmonization ─────────────────────────────────────────────────
            % Shared target spec for DeePhys and external data.
            % buildFeatureMatrix reads these to ensure DeePhys and external units
            % are projected into the same feature space.
            defaultParams.Harmonization.WaveformTargetSamplingRate = 120000;  % Hz
            defaultParams.Harmonization.WaveformPreTrough          = 1.0;     % ms before trough
            defaultParams.Harmonization.WaveformPostTrough         = 2.0;     % ms after trough
            % WaveformEdgeMode: "edge" (default) pads short waveforms with their edge values so
            %   every unit (and every external dataset) gets the same fixed window; "trim"
            %   shrinks the window to the shortest span of the units at hand, so training and
            %   external feature vectors could differ in length; "zero" pads with zeros.
            defaultParams.Harmonization.WaveformEdgeMode           = "edge";
            defaultParams.Harmonization.ACGBinSize                 = 0.0005;  % s (0.5 ms bins)
            defaultParams.Harmonization.ACGLag                     = 0.1;     % s (+-100 ms -> 401 bins)
            % ACGSource: "ACG" (the unit's baseline segment only; default), "FullACG"
            %   (stored parent ACG) or "Concatenated" (all segments joined). With drug-response
            %   labels, "ACG" keeps the features drug-free: the other two include the drug
            %   periods, so responsive units' ACGs encode the response used to label them.
            %   "FullACG" is appropriate when the segments share one state (no perturbation),
            %   where it only adds spikes.
            defaultParams.Harmonization.ACGSource                  = "ACG";
            % ACGSmoothing: Gaussian width (ms) applied to the ACGs before feature building,
            %   or "auto" (default): chosen without labels by split-half prediction
            %   (selectACGSmoothing). External units get the same width. 0 = none.
            defaultParams.Harmonization.ACGSmoothing               = "auto";
            % ACGNormalization: "peak" (default; divide by the maximum) or "baseline" (divide
            %   by the mean of the outer fifth of the lags: 1 = uncorrelated level, the shape
            %   no longer scales with firing rate / burstiness).
            defaultParams.Harmonization.ACGNormalization           = "peak";
        end

    end

    % =====================================================================
    % Diagnostic methods (callable independently or auto-triggered)
    % =====================================================================
    methods

        % Declared here; implemented in diagnosticResponsiveUnits.m
        diagnosticResponsiveUnits(ctc)

        % Declared here; implemented in diagnosticTrainLabels.m
        diagnosticTrainLabels(ctc)

        % Declared here; implemented in diagnosticClassification.m
        diagnosticClassification(ctc)

        % Declared here; implemented in diagnosticOptimization.m
        diagnosticOptimization(ctc, results)

        % Declared here; implemented in evaluatePipeline.m
        cv_results = evaluatePipeline(ctc, opts)

    end

    % =====================================================================
    % Private helpers
    % =====================================================================
    methods (Access = private)

        % Declared here; implemented in buildNormalizedFeatures.m
        buildNormalizedFeatures(ctc)

        % Declared here; implemented in saveDiagnosticFigure.m
        saveDiagnosticFigure(ctc, fig, stage_name)

        function [unique_ud, all_to_unique, unique_to_all] = uniqueUnitMap(ctc)
        % UNIQUEUNITMAP  Return one representative UnitData per unique UnitID.
        %
        % For cultures with multiple recordings (e.g. dose-response), the same
        % physical unit appears once per recording in UnitDataArray. Classification
        % operates on unique units; results are then broadcast back to all rows.
        %
        % The representative recording is the one whose GroupingVar value matches
        % Parameters.UMAP.GroupingValues (the baseline). Falls back to the first
        % occurrence if no baseline recording is found.
        %
        % OUTPUTS:
        %   unique_ud      - UnitData array, one entry per unique UnitID
        %   all_to_unique  - (1 x numel(UnitDataArray)) index into unique_ud
        %   unique_to_all  - (1 x numel(unique_ud)) cell of indices into UnitDataArray
            ud          = ctc.UnitDataArray;
            uid_strings = string({ud.UnitID})';
            [unique_uids, ~, all_to_unique] = unique(uid_strings, 'stable');
            n_unique    = numel(unique_uids);

            meta         = ctc.FeatureStore.MetadataTable;
            group_field  = string(ctc.Parameters.UMAP.GroupingVar);
            baseline_gv  = ctc.Parameters.UMAP.GroupingValues;
            has_group    = ismember(group_field, meta.Properties.VariableNames);

            unique_to_all   = cell(1, n_unique);
            best_idx        = zeros(1, n_unique);

            for u = 1:n_unique
                rows = find(all_to_unique == u);
                unique_to_all{u} = rows';
                best_idx(u) = rows(1);  % default: first occurrence

                if ~has_group; continue; end
                for ri = 1:numel(rows)
                    rec_row = meta(meta.RecordingID == string(ud(rows(ri)).RecordingID), :);
                    if isempty(rec_row); continue; end
                    gv = rec_row.(group_field)(1);
                    if isnumeric(gv) && isnumeric(baseline_gv) && gv == baseline_gv
                        best_idx(u) = rows(ri);
                        break
                    elseif (ischar(gv) || isstring(gv)) && string(gv) == string(baseline_gv)
                        best_idx(u) = rows(ri);
                        break
                    end
                end
            end

            unique_ud     = ud(best_idx);
            all_to_unique = all_to_unique';  % (1 x N)
        end

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
                sr = [unit_data.SamplingRate];        % original sampling rate, per unit
            end

            % ACG: compute from spike times at the specified bin size / lag
            n_bins = round(ph.ACGLag / ph.ACGBinSize) * 2 + 1;
            acg    = zeros(n_bins, N);
            for u = 1:N
                acg(:, u) = CellTypeClassifier.computeACG( ...
                    unit_data(u).SpikeTimes, ph.ACGBinSize, ph.ACGLag);
            end
        end

        function acg = computeParentACGs(ctc, unit_data, ph)
        % COMPUTEPARENTACGS  ACGs from concatenated spike trains across all recording segments.
        %
        % When OriginalUnitDataArray holds multiple recordings per unit (dose-response),
        % concatenates spike trains from all segments of each unit (with inter-segment
        % gaps to prevent cross-segment ISIs from contaminating the ACG).
            orig_ud         = ctc.OriginalUnitDataArray;
            all_uid_strings = string({orig_ud.UnitID});
            n_units = numel(unit_data);
            n_bins  = round(ph.ACGLag / ph.ACGBinSize) * 2 + 1;
            acg     = zeros(n_bins, n_units);
            gap     = 2 * ph.ACGLag;   % prevents cross-segment ISIs within ACG window

            for u = 1:n_units
                uid         = string(unit_data(u).UnitID);
                segment_idx = find(all_uid_strings == uid);

                if numel(segment_idx) <= 1
                    acg(:, u) = CellTypeClassifier.computeACG( ...
                        unit_data(u).SpikeTimes, ph.ACGBinSize, ph.ACGLag);
                else
                    combined_spikes = [];
                    offset = 0;
                    for si = 1:numel(segment_idx)
                        st = orig_ud(segment_idx(si)).SpikeTimes(:);
                        if ~isempty(st)
                            combined_spikes = [combined_spikes; st + offset]; %#ok<AGROW>
                        end
                        offset = offset + orig_ud(segment_idx(si)).RecordingDuration + gap;
                    end
                    acg(:, u) = CellTypeClassifier.computeACG( ...
                        combined_spikes, ph.ACGBinSize, ph.ACGLag);
                end
            end
        end

    end

    methods (Static)

        function binned = binnedSpikeMatrix(unit_data, time_window, bin_size)
        % BINNEDSPIKEMATRIX  (N_units x N_bins) spike count matrix from UnitData.SpikeTimes.
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
        end

        function W = knnGraph(X, k)
        % KNNGRAPH  Deterministic UMAP-style fuzzy kNN graph with exact neighbours.
        %   W = CellTypeClassifier.knnGraph(X, k)   (see private/buildKNNGraph)
            W = buildKNNGraph(X, k);
        end

        function [W, info] = wnnGraph(X, is_acg, k)
        % WNNGRAPH  Weighted-nearest-neighbour graph of the ACG and waveform columns.
        %   [W, info] = CellTypeClassifier.wnnGraph(X, is_acg, k)   (see private/buildWNNGraph)
            [W, info] = buildWNNGraph(X, is_acg, k);
        end

        function [sigma_ms, info] = acgSmoothingWidth(spike_trains, acg_used, bin_size, lag)
        % ACGSMOOTHINGWIDTH  Label-free ACG smoothing width (see private/selectACGSmoothing).
            [sigma_ms, info] = selectACGSmoothing(spike_trains, acg_used, bin_size, lag);
        end

        function [t, info] = externalCutoff(p, min_dbic)
        % EXTERNALCUTOFF  Label-free threshold on external P(inhibitory) (see private/inferExternalCutoff).
            if nargin < 2, min_dbic = 10; end
            [t, info] = inferExternalCutoff(p, min_dbic);
        end

        function acg = smoothACGs(acg, sigma_ms, bin_size)
        % SMOOTHACGS  Gaussian smoothing of ACG columns (centre bin kept; see private/smoothACG).
            acg = smoothACG(acg, sigma_ms / 1000 / bin_size);
        end

        function T = responseTable(ud_all, meta, culture_keys, group_field, baseline_gv, cultures, p)
        % RESPONSETABLE  Drift-corrected baseline vs high-dose response per unit
        %   (see private/calibratedResponseTable).
            T = calibratedResponseTable(ud_all, meta, culture_keys, group_field, baseline_gv, cultures, p);
        end

        function [inc, dec, info] = calibrateResponses(T, is_ctrl, q, alpha, target_fdr)
        % CALIBRATERESPONSES  Responder calls with control-culture thresholds (see private/calibrateCalls).
            if nargin < 5, target_fdr = NaN; end
            [inc, dec, info] = calibrateCalls(T, is_ctrl, q, alpha, target_fdr);
        end

        function [T, inh_comm, exc_comm] = communityEnrichment(M, eligible, responsive, fdr)
        % COMMUNITYENRICHMENT  Hypergeometric responder enrichment/depletion per community.
        %   [T, inh, exc] = CellTypeClassifier.communityEnrichment(M, eligible, responsive, fdr)
        %   (see private/communityEnrichmentTest)
            [T, inh_comm, exc_comm] = communityEnrichmentTest(M, eligible, responsive, fdr);
        end

        function [inh_seeds, exc_seeds] = enrichmentSeeds(G, responsive, eligible, p_comm, rng_seed)
        % ENRICHMENTSEEDS  Training seeds of the enrichment test on graph G, all units in the
        %   training subset (see private/selectEnrichmentSeeds). Returns global unit indices.
            if nargin < 5, rng_seed = 1; end
            responsive = logical(responsive(:)');
            subset = 1:size(G, 1);
            [~, det] = evalc('selectEnrichmentSeeds(G, find(responsive), subset, responsive, logical(eligible(:)''), p_comm, rng_seed)');
            inh_seeds = subset(det.clean_resp_local);
            exc_seeds = subset(det.counterexample_idx_local);
        end

        function [labels, kept, expected_error] = selectConfident(p_inh, n_spikes, min_p, min_spikes)
        % SELECTCONFIDENT  Labels of units with P(inhibitory) >= min_p or <= 1 - min_p and at
        %   least min_spikes spikes; NaN otherwise (see confidentLabels). expected_error is
        %   mean(1 - max(P, 1 - P)) over the kept units, the share of wrong labels if P is
        %   calibrated.
            p_inh = p_inh(:)'; n_spikes = n_spikes(:)';
            kept = (p_inh >= min_p | p_inh <= 1 - min_p) & n_spikes >= min_spikes & ~isnan(p_inh);
            labels = nan(size(p_inh));
            labels(kept) = 1 + double(p_inh(kept) > 0.5);
            expected_error = mean(min(p_inh(kept), 1 - p_inh(kept)));
            if ~any(kept), expected_error = NaN; end
        end

        function rules = parseResponseRules(spec)
        % PARSERESPONSERULES  Bootstrap.ResponseRules as a struct array (Filter = {field, value},
        %   Direction). spec: n x 3 cell {field, value, direction} or a struct array with fields
        %   Filter (or Field + Value) and Direction.
            if iscell(spec)
                assert(size(spec, 2) == 3, 'CellTypeClassifier:badResponseRules', ...
                    'ResponseRules must have one row {field, value, direction} per rule.');
                rules = struct('Filter', cell(1, size(spec, 1)), 'Direction', "");
                for r = 1:size(spec, 1)
                    rules(r).Filter = spec(r, 1:2);
                    rules(r).Direction = lower(string(spec{r, 3}));
                end
            else
                rules = struct('Filter', cell(1, numel(spec)), 'Direction', "");
                for r = 1:numel(spec)
                    if isfield(spec, 'Filter'), f = spec(r).Filter; else, f = {spec(r).Field, spec(r).Value}; end
                    rules(r).Filter = f;
                    rules(r).Direction = lower(string(spec(r).Direction));
                end
            end
            for r = 1:numel(rules)
                assert(ismember(rules(r).Direction, ["increase", "decrease", "both"]), ...
                    'CellTypeClassifier:badResponseRules', ...
                    'Rule %d: direction must be "increase", "decrease" or "both".', r);
            end
        end

        function st = thinSpikeTrain(st, n, rs)
        % THINSPIKETRAIN  Random subset of n spikes (sorted); all spikes if n >= numel(st).
            if nargin < 3, rs = RandStream.getGlobalStream; end
            if n < numel(st)
                st = sort(st(randperm(rs, numel(st), n)));
            end
        end

        function acg = computeACG(spike_times, bin_size, lag, sampling_rate)
        % COMPUTEACG  Autocorrelogram for a single unit from spike times.
        %
        %   acg = CellTypeClassifier.computeACG(spike_times)
        %   acg = CellTypeClassifier.computeACG(spike_times, bin_size, lag)
        %   acg = CellTypeClassifier.computeACG(spike_times, bin_size, lag, sampling_rate)
        %
        %   spike_times   - (N x 1) spike times in seconds
        %   bin_size      - bin width in seconds (default: 0.0005)
        %   lag           - half-width of the ACG window in seconds (default: 0.1)
        %   sampling_rate - recording sampling rate in Hz (default: 20000, Maxwell HDMEA)
        %
        %   Returns acg as (N_bins x 1).  Use computeACGs for batches.
            arguments
                spike_times   (:,1) double
                bin_size      (1,1) double = 0.0005
                lag           (1,1) double = 0.1
                sampling_rate (1,1) double = 20000
            end
            n_bins = round(lag / bin_size) * 2 + 1;
            acg    = zeros(n_bins, 1);
            if numel(spike_times) > 1
                [acg_u, ~] = CCG(spike_times, ones(numel(spike_times), 1), ...
                    'binSize', bin_size, 'duration', lag * 2, 'Fs', 1 / sampling_rate);
                if numel(acg_u) == n_bins
                    acg = acg_u;
                end
            end
        end

        function acgs = computeACGs(spike_times_cell, bin_size, lag, sampling_rate)
        % COMPUTEACGS  Autocorrelograms for multiple units from spike times.
        %
        %   acgs = CellTypeClassifier.computeACGs(spike_times_cell)
        %   acgs = CellTypeClassifier.computeACGs(spike_times_cell, bin_size, lag)
        %   acgs = CellTypeClassifier.computeACGs(spike_times_cell, bin_size, lag, sampling_rate)
        %
        %   spike_times_cell - (1 x N) cell array of spike-time vectors (seconds)
        %   bin_size         - bin width in seconds (default: 0.0005)
        %   lag              - half-width of the ACG window in seconds (default: 0.1)
        %   sampling_rate    - recording sampling rate in Hz (default: 20000, Maxwell HDMEA)
        %
        %   Returns acgs as (N_bins x N), ready for classifyExternalUnits.
            arguments
                spike_times_cell (1,:) cell
                bin_size         (1,1) double = 0.0005
                lag              (1,1) double = 0.1
                sampling_rate    (1,1) double = 20000
            end
            N    = numel(spike_times_cell);
            acgs = zeros(round(lag / bin_size) * 2 + 1, N);
            for u = 1:N
                acgs(:, u) = CellTypeClassifier.computeACG( ...
                    spike_times_cell{u}(:), bin_size, lag, sampling_rate);
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
            assert(numel(y_pred) == numel(y_true), 'CellTypeClassifier:lengthMismatch', ...
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

        function [bin_size, lag] = storedACGParams(unit_data, n_bins)
        % STOREDACGPARAMS  Bin size and lag of stored parent ACGs ([] if unknown).
        %   RecordingProcessor.computeParentFeatures records them on each unit; older
        %   processors without that provenance are assumed to use its defaults when the
        %   bin count matches them (0.5 ms bins, +-100 ms -> 401 bins).
            bin_size = []; lag = [];
            if ~isempty(unit_data) && isprop(unit_data(1), 'FullACGBinSize') && ...
                    ~isempty(unit_data(1).FullACGBinSize) && ~isempty(unit_data(1).FullACGLag)
                bin_size = unit_data(1).FullACGBinSize;
                lag      = unit_data(1).FullACGLag;
            elseif n_bins == 401
                bin_size = 0.0005;
                lag      = 0.1;
            end
        end

        function out = adaptACG(acg, from_bs, from_lag, to_bs, to_lag)
        % ADAPTACG  Crop (and re-bin by an integer factor) bin-centred ACGs to a narrower,
        %   coarser window. Returns [] when the target cannot be derived from the source.
            out = [];
            if isempty(from_bs) || to_lag > from_lag + 1e-12
                return
            end
            k = to_bs / from_bs;
            if k < 1 - 1e-9 || abs(k - round(k)) > 1e-9
                return
            end
            k = round(k);
            n_from = size(acg, 1);
            if n_from ~= round(2 * from_lag / from_bs) + 1
                return
            end
            if k > 1
                % A coarse bin centred on a fine-bin centre spans k fine bins; for even k
                % the two outermost fine bins are split in half.
                if mod(k, 2) == 0
                    kernel = [0.5, ones(1, k - 1), 0.5];
                else
                    kernel = ones(1, k);
                end
                acg = conv2(acg, kernel(:), 'same');
            end
            centre = (n_from + 1) / 2;
            m = round(to_lag / to_bs);
            out = acg(centre + (-m:m) * k, :);
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
