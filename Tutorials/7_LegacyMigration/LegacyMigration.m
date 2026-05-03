%% Tutorial 7 — Legacy Migration
%
% Converts old MEArecording / RecordingGroup objects to the new API.
% After migration, all downstream analyses (Tutorials 2–6) work without change.
%
% Migration paths:
%   A. Single recording  : RecordingProcessor.fromLegacyMat(path)
%   B. RecordingGroup    : Experiment.fromLegacyGroup(rg)
%   C. Raw MEArecording  : FeatureStore.fromLegacyRecordings(mearec_array)
%   D. Per-unit          : UnitData.fromLegacyUnit(unit_obj, recording_id)
%
% Prerequisites:
%   - Old .mat files with MEArecording or RecordingGroup objects
%   - DeePhys on the MATLAB path

%% ============================================================
%% A  Single recording migration
%% ============================================================

%% A1  Convert MEArecording.mat to RecordingProcessor

old_mat  = '/path/to/old/MEArecording.mat';
save_dir = '/path/to/new';

proc = RecordingProcessor.fromLegacyMat(old_mat);

fprintf('Migrated %d units from %s\n', numel(proc.Units), old_mat);
fprintf('Status: QC=%s UnitFeatures=%s NetworkFeatures=%s\n', ...
    proc.Status.QC, proc.Status.UnitFeatures, proc.Status.NetworkFeatures);

%% A2  Inspect migrated features

if ~isempty(proc.UnitFeatureTable)
    fprintf('UnitFeatureTable: %d units × %d features\n', ...
        height(proc.UnitFeatureTable), width(proc.UnitFeatureTable));
end

if ~isempty(proc.NetworkFeatureTable)
    disp(proc.NetworkFeatureTable);
end

%% A3  Save in new format

proc_file = fullfile(save_dir, 'RecordingProcessor.mat');
proc.save(proc_file);

% Verify round-trip
proc2 = RecordingProcessor.load(proc_file);
assert(numel(proc2.Units) == numel(proc.Units), 'Unit count mismatch after save/load');

%% A4  Batch migration

old_mats = {'/path/to/rec1.mat', '/path/to/rec2.mat', '/path/to/rec3.mat'};
procs(numel(old_mats)) = RecordingProcessor();
for i = 1:numel(old_mats)
    try
        procs(i) = RecordingProcessor.fromLegacyMat(string(old_mats{i}));
    catch ME
        warning('Migration failed for %s: %s', old_mats{i}, ME.message);
    end
end
fprintf('Successfully migrated %d / %d recordings\n', ...
    sum(arrayfun(@(p) ~isempty(p.Units), procs)), numel(procs));

%% ============================================================
%% B  RecordingGroup migration
%% ============================================================

%% B1  Load old RecordingGroup

load('/path/to/RecordingGroup.mat', 'rg');   % variable name may differ

%% B2  Build Experiment from RecordingGroup

exp = Experiment.fromLegacyGroup(rg);

fprintf('FeatureStore: %d units, %d recordings\n', ...
    height(exp.FeatureStore.UnitTable), height(exp.FeatureStore.RecordingTable));

%% B3  Run analyses on migrated Experiment (same as Tutorial 3)

result = exp.classify('Unit', 'Mutation', struct('Algorithm', 'rf', 'KFold', 5));
summary = ClassificationResult.summarizeFolds(result);
fprintf('Accuracy: %.2f +/- %.2f\n', summary.mean_accuracy, summary.std_accuracy);

%% ============================================================
%% C  Direct FeatureStore migration from MEArecording array
%% ============================================================

%% C1  Load old objects

% old_recs is an array of MEArecording objects (already in memory or loaded)
%   old_recs = [mearec1, mearec2, mearec3];

% For a RecordingGroup:
%   old_recs = rg.Recordings;

%% C2  Build FeatureStore directly

% fs = FeatureStore.fromLegacyRecordings(old_recs);
%
% Optionally select feature groups:
%   fs = FeatureStore.fromLegacyRecordings(old_recs, "all", "all");

% Once migrated, save and use as usual:
%   fs.save('/path/to/FeatureStore.mat');

%% ============================================================
%% D  Per-unit migration
%% ============================================================

%% D1  Convert a Unit object to UnitData

% unit_obj is an old Unit handle object
% rec_id   = paramHash(struct('InputPath', char(unit_obj.MEArecording.Metadata.InputPath)));
%
% ud = UnitData.fromLegacyUnit(unit_obj, rec_id);
%
% Bulk conversion from Unit array:
%   ud_array = UnitData.fromLegacyUnitArray(unit_array, rec_id);

%% ============================================================
%% E  Verification
%% ============================================================

%% E1  Compare old and new feature values
%
% Features should be numerically identical if the same parameters were used.

if exist('proc', 'var') && exist('old_mat', 'var')
    try
        data_old = load(old_mat);
        fns = fieldnames(data_old);
        old_obj = [];
        for i = 1:numel(fns)
            if isa(data_old.(fns{i}), 'MEArecording')
                old_obj = data_old.(fns{i});
                break
            end
        end

        if ~isempty(old_obj) && ~isempty(proc.UnitFeatureTable) && ~isempty(old_obj.Units)
            fprintf('Old units: %d, New units: %d\n', numel(old_obj.Units), numel(proc.Units));
        end
    catch
        warning('Could not load old object for verification.');
    end
end

%% E2  Smoke-test: FeatureStore extraction works

if exist('exp', 'var')
    [X, ~] = exp.FeatureStore.unitMatrix('all');
    fprintf('Feature matrix: %d × %d (%.1f%% NaN)\n', ...
        height(X), width(X), 100 * mean(isnan(table2array(X)), 'all'));
end
