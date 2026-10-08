function generateTrainLabels(ctc, opts)
% GENERATETRAINLABELS  Build training labels via outlier detection and feature-space operations.
%
% Pipeline:
%   1. Global normalization + optional feature selection  (computeGlobalNormalization)
%   2. Unsupervised UMAP on all unique units              (runUnsupervisedUMAP)
%   3. Outlier detection and CE selection                 (selectExplicitCE |
%                                                          selectEnrichmentSeeds (default) |
%                                                          detectCommunityOutliers |
%                                                          selectIforestCE | selectNoCE)
%   4. Assemble training labels                           (assembleTrainLabels)
%
% Two label-source strategies:
%   Drug response  — responsive units = inhibitory ground truth; CEs inferred.
%   Metadata       — ctc.CounterexampleUnitIdx provides explicit excitatory labels. With
%                    OutlierDetection.Method "enrichment" (default) the seeds come from the same
%                    enrichment test (labelled units = eligible pool); other methods use the
%                    iforest cleaning of selectExplicitCE.
%
% Requires ctc.ResponsiveUnitIdx (run identifyResponsiveUnits first).
% Sets: ctc.TrainLabels, ctc.NormalizationParams, ctc.Reduction.Unsupervised

arguments
    ctc CellTypeClassifier
    opts.ReuseExistingUMAP (1,1) logical = true
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before generateTrainLabels()');

p_umap  = ctc.Parameters.UMAP;
p_outlr = ctc.Parameters.OutlierDetection;
p_train = ctc.Parameters.TrainLabels;

rng(ctc.Parameters.RNGSeed, 'twister');
ctc.clearCache();
ctc.buildNormalizedFeatures();
nf = ctc.NormalizedFeatures;

% -- Stage 1: global normalization + feature selection ------------------------
[X_all, norm_params] = computeGlobalNormalization(nf, p_umap, ctc.Parameters.Harmonization);
norm_params.acg_sigma_ms = nf.acg_sigma_ms;   % applied to external ACGs too
ctc.NormalizationParams = norm_params;

% -- Stage 2: unsupervised UMAP -----------------------------------------------
[reduction, n_neighbors] = ctc.runUnsupervisedUMAP(X_all, p_umap, opts.ReuseExistingUMAP);

% -- Map responsive units to training-subset local indices --------------------
unit_ids_table    = string(ctc.FeatureStore.UnitTable.UnitID);
unique_ud_ids     = string({nf.unique_ud.UnitID});
resp_uids         = unique(unit_ids_table(ctc.ResponsiveUnitIdx));
resp_unique       = ismember(unique_ud_ids, resp_uids);
subset_responsive = resp_unique(nf.subset_mask);
responsive_local  = find(subset_responsive);
subset_global_idx = find(nf.subset_mask);

if isempty(responsive_local)
    n_resp_total = sum(ctc.ResponsiveUnitIdx);
    error('CellTypeClassifier:noResponsiveUnits', ...
        ['No responsive units found in the training culture subset.\n' ...
         '  Total responsive units across all cultures: %d\n' ...
         '  Training subset size: %d unique units\n' ...
         'Possible causes:\n' ...
         '  1. Parameters.UMAP.TrainingCultureIdx excludes cultures with responsive units.\n' ...
         '  2. Parameters.UMAP.GroupingValues does not match the baseline dose in the data.\n' ...
         '  3. identifyResponsiveUnits() was not run, or produced 0 responsive units.\n' ...
         'Diagnostic: check sum(ctc.ResponsiveUnitIdx) and ctc.Parameters.UMAP.TrainingCultureIdx.'], ...
        n_resp_total, sum(nf.subset_mask));
end

X_feat          = X_all(nf.subset_mask, :);
has_explicit_ce = ~isempty(ctc.CounterexampleUnitIdx) && any(ctc.CounterexampleUnitIdx);

% -- Stage 3: outlier detection and CE selection ------------------------------
if has_explicit_ce && lower(string(p_outlr.Method)) == "enrichment"
    % Both classes labelled (e.g. pure cultures): the same enrichment test, with the labelled
    % units as the eligible pool. Seeds are labelled units whose community is enriched for
    % their own class, so contaminants of the other type in a "pure" culture are left out.
    labelled_fs = ctc.ResponsiveUnitIdx(:) | ctc.CounterexampleUnitIdx(:);
    eligible_global = false(1, nf.n_unique);
    eligible_global(:) = labelled_fs(nf.unique_to_fs);
    det = selectEnrichmentSeeds(ctc.graphMatrix(), responsive_local, subset_global_idx, ...
        resp_unique, eligible_global, ctc.Parameters.Community, ctc.Parameters.RNGSeed);
elseif has_explicit_ce
    ce_uid_table = unique(unit_ids_table(ctc.CounterexampleUnitIdx));
    ce_subset    = ismember(unique_ud_ids(nf.subset_mask), ce_uid_table);
    ce_local_all = find(ce_subset);
    det = selectExplicitCE(X_feat, responsive_local, subset_responsive, ...
        subset_global_idx, ce_local_all, reduction, p_outlr);
else
    odm = lower(string(p_outlr.Method));

    comm_diagnostics = [];
    if odm == "enrichment"
        % Eligible units: those in cultures where responders occur (e.g. hM3Dq cultures),
        % so the enrichment test is not confounded by culture condition.
        ut      = ctc.FeatureStore.UnitTable;
        cult    = FeatureStore.getCultureIDsForUnits(ut.RecordingID, ctc.FeatureStore.MetadataTable, ...
                      ctc.Parameters.CultureKeys);
        resp_cultures = unique(cult(ctc.ResponsiveUnitIdx));
        elig_fs = ismember(cult, resp_cultures);
        eligible_global = false(1, nf.n_unique);
        eligible_global(:) = elig_fs(nf.unique_to_fs);
        det = selectEnrichmentSeeds(ctc.graphMatrix(), responsive_local, subset_global_idx, ...
            resp_unique, eligible_global, ctc.Parameters.Community, ctc.Parameters.RNGSeed);
    elseif odm == "community"
        det = detectCommunityOutliers(X_feat, responsive_local, subset_responsive, ...
            subset_global_idx, reduction, size(X_all, 1), n_neighbors, resp_unique, ...
            ctc.Parameters.Community, ctc.graphMatrix(), ctc.Parameters.RNGSeed);
        if ~det.success
            comm_diagnostics = det;   % preserve Louvain results even on fallback
            odm = "iforest";
        end
    end

    if odm == "iforest"
        det = selectIforestCE(X_feat, responsive_local, subset_responsive, ...
            subset_global_idx, reduction, p_outlr);
        if ~isempty(comm_diagnostics)
            det.community_ids = comm_diagnostics.community_ids;
            det.inh_comm_ids  = comm_diagnostics.inh_comm_ids;
            det.Q_modularity  = comm_diagnostics.Q_modularity;
        end
    elseif odm == "none"
        det = selectNoCE(X_feat, responsive_local, subset_responsive, p_outlr);
    end
end

det.has_explicit_ce = has_explicit_ce;

% -- Stage 4: assemble training labels ----------------------------------------
ctc.TrainLabels = assembleTrainLabels(det, subset_global_idx, nf, p_train, ...
    unit_ids_table, ctc.ResponsiveStrength);

resp_label    = p_train.ResponsiveClassLabel;
counter_label = 3 - resp_label;
fprintf('Training set: %i %s, %i %s candidates\n', ...
    sum(ctc.TrainLabels.sorted_y_train == counter_label), localLabelName(counter_label), ...
    sum(ctc.TrainLabels.sorted_y_train == resp_label),    localLabelName(resp_label));

if ctc.Parameters.Diagnostics.Enable
    ctc.diagnosticTrainLabels();
end
end

function name = localLabelName(label)
    if label == 1; name = 'excitatory';
    elseif label == 2; name = 'inhibitory';
    else; name = sprintf('class_%d', label);
    end
end
