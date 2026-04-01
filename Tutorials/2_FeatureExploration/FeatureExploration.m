%% Tutorial 2 — Feature Exploration
%
% Demonstrates how to inspect and extract feature matrices from a FeatureStore
% for plotting, quality checking, and preparing ML inputs.
%
% Prerequisites:
%   - Tutorial 1 completed (FeatureStore saved to disk)
%   - DeePhys on the MATLAB path

%% 1  Load a saved FeatureStore

fs_file = '/path/to/FeatureStore.mat';
fs = FeatureStore.load(fs_file);

fprintf('Units     : %d\n', height(fs.UnitTable));
fprintf('Recordings: %d\n', height(fs.RecordingTable));

%% 2  Inspect table columns

% All column names in UnitTable
all_cols = string(fs.UnitTable.Properties.VariableNames)';
disp(all_cols);

% Separate identity / metadata / feature columns
id_cols   = ["UnitID", "RecordingID"];
meta_cols = string(fs.MetadataTable.Properties.VariableNames(2:end))';  % skip RecordingID
fprintf('\nMetadata columns:\n');  disp(meta_cols);

feat_cols = all_cols(~ismember(all_cols, [id_cols; meta_cols]));
fprintf('Feature columns: %d total\n', numel(feat_cols));

% Parent feature columns (from parent recording, stored with prefix)
parent_feat_cols = feat_cols(startsWith(feat_cols, 'Parent_'));
child_feat_cols  = feat_cols(~startsWith(feat_cols, 'Parent_'));
fprintf('  Child features : %d\n', numel(child_feat_cols));
fprintf('  Parent features: %d\n', numel(parent_feat_cols));

%% 3  Extract a full unit feature matrix

% "all" returns only child feature columns (Parent_* excluded)
[X_all, unit_ids] = fs.unitMatrix('all');
fprintf('unitMatrix("all"): %d units × %d features\n', height(X_all), width(X_all));

%% 4  Select specific feature groups

[X_act, ~]  = fs.unitMatrix('ActivityFeatures');
[X_wf, ~]   = fs.unitMatrix('WaveformFeatures');
[X_acg, ~]  = fs.unitMatrix('ACG');
[X_multi, ~] = fs.unitMatrix(["ActivityFeatures", "WaveformFeatures"]);

fprintf('ActivityFeatures  : %d cols\n', width(X_act));
fprintf('WaveformFeatures  : %d cols\n', width(X_wf));
fprintf('ACG               : %d cols\n', width(X_acg));
fprintf('Activity+Waveform : %d cols\n', width(X_multi));

%% 5  Use parent features instead of child features
%
% Parent ACGs (from the full concatenated recording) generally give better
% estimates because more spikes → smoother autocorrelogram.
% Output column names are always child names (ACG1, not Parent_ACG1) so
% feature matrices from datasets with/without parent features are compatible.

% Prefer parent ACG, fall back to child ACG if Parent_ACG* not present
[X_parent_acg, ~] = fs.unitMatrix('all', 'ACG');
fprintf('unitMatrix with parent ACG: %d × %d\n', height(X_parent_acg), width(X_parent_acg));

% Use all available parent features
[X_parent_all, ~] = fs.unitMatrix('all', 'all');

%% 6  Recording-level feature matrix

[X_rec, rec_ids] = fs.recordingMatrix('all');
fprintf('recordingMatrix("all"): %d recordings × %d features\n', height(X_rec), width(X_rec));

%% 7  Culture-level feature matrix
%
% Groups recordings by culture identity (ChipID + PlatingDate), selects
% recordings at each specified DIV, and returns a wide table with
% feature columns suffixed by grouping value (e.g., FiringRate_7, FiringRate_14).

identity_keys   = ["ChipID", "PlatingDate"];
grouping_var    = "DIV";
grouping_values = [7, 14, 21, 28];

[X_culture, culture_ids] = fs.cultureMatrix(identity_keys, grouping_var, grouping_values);
fprintf('cultureMatrix: %d cultures × %d cols\n', height(X_culture), width(X_culture));
disp(culture_ids);

% With normalization: divide each time point by the DIV 7 baseline
[X_norm, ~] = fs.cultureMatrix(identity_keys, grouping_var, grouping_values, 'baseline');

%% 8  Subset by metadata

% Keep only WT recordings
fs_wt = fs.subsetByMetadata('Mutation', 'WT');
fprintf('WT units: %d / %d\n', height(fs_wt.UnitTable), height(fs.UnitTable));

% Keep recordings at DIV 21
fs_d21 = fs.subsetByMetadata('DIV', 21);

% Keep multiple values
fs_ko_wt = fs.subsetByMetadata('Mutation', {'WT', 'KO'});

%% 9  Plot waveform overlays

wf_cols = startsWith(string(X_wf.Properties.VariableNames), 'Waveform');
if any(wf_cols)
    wf_mat = X_wf.Variables;
    figure;
    plot(wf_mat(1:min(50, height(wf_mat)), :)', 'Color', [0 0 0 0.3]);
    xlabel('Sample'); ylabel('Amplitude (norm.)');
    title('Reference waveforms (first 50 units)');
end

%% 10  Plot ACG distributions

if width(X_acg) > 0
    acg_mat = X_acg.Variables;
    n_bins  = width(X_acg);
    lag_ms  = linspace(-100, 100, n_bins);
    figure;
    imagesc(lag_ms, 1:min(100, height(acg_mat)), acg_mat(1:min(100, end), :));
    xlabel('Lag (ms)'); ylabel('Unit index');
    title('ACGs (first 100 units)');
    colorbar;
end

%% 11  Feature distributions by metadata group

if ismember('Mutation', string(fs.UnitTable.Properties.VariableNames))
    mutations = unique(string(fs.UnitTable.Mutation));
    if width(X_act) > 0
        fr_col = 'MeanFiringRate';
        if ismember(fr_col, X_act.Properties.VariableNames)
            figure; hold on;
            for m = 1:numel(mutations)
                mask = string(fs.UnitTable.Mutation) == mutations(m);
                histogram(X_act.(fr_col)(mask), 'Normalization', 'pdf', ...
                    'DisplayName', mutations(m));
            end
            xlabel('Mean firing rate (Hz)'); ylabel('Density');
            legend; title('Firing rate by genotype');
        end
    end
end

%% 12  Missing values

X_mat = table2array(fs.unitMatrix('all'));
nan_frac_col = mean(isnan(X_mat), 1);
nan_frac_row = mean(isnan(X_mat), 2);
fprintf('Features with >10%% NaN: %d / %d\n', sum(nan_frac_col > 0.1), numel(nan_frac_col));
fprintf('Units with any NaN    : %d / %d\n', sum(nan_frac_row > 0),   numel(nan_frac_row));
