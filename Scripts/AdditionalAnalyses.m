%% If you have to load the files first
failed_list = [];
parfor iRec = 1:length(sorting_path_list)
    mearec = fullfile(sorting_path_list(iRec), "MEArecording.mat");
    if exist(mearec,'file')
        try
            opened = load(mearec,"obj");
            mearec = opened.obj;
            if length([mearec.Units]) >= 20
                mearec.Parameters.Save.Overwrite = true; %We need to save the objects to use parfor loops // properties are not changed
                mearec.inferGraphFeatures(); % Select analyses that you want to run
                mearec.UnitFeatures = mearec.aggregateSingleCellFeatures(mearec.Units);
                mearec.saveObject(); % Save
            end
        catch
            failed_list = [failed_list  mearec];
        end
    end
end

%% If you already have them in your workspace
parfor iRec = 1:length(rec_array)
%     if ~isfield(rec_array(iRec).NetworkFeatures,'GraphFeatures') || isempty(rec_array(iRec).NetworkFeatures.GraphFeatures)
        rec_array(iRec).generateUnits();
        rec_array(iRec).saveObject();
%     end
end


%%

for r = 1:length(rec_array)
   metadata = rec_array(r).Metadata;
   metadata.PathPartIndex = [11,13,14];
   rec_array(r).parseMetadata(metadata);
end
%%
iMissing = arrayfun(@(x) isempty(x.Catch22),rec_group.Units);
%%
missing_units = rec_group.Units(iMissing);
for m = 1:length(missing_units)
   catch_22_table = missing_units(m).MEArecording.run_catch_22(missing_units(m).SpikeTimes);
   catch_22_table.Properties.VariableNames = "SC_" + string(catch_22_table.Properties.VariableNames);
   missing_units(m).Catch22 = catch_22_table;
end

%%
parfor r = 189:length(rec_array)
    rec_array(r).Parameters.Save.Overwrite = true; 
    rec_array(r).UnitFeatures = rec_array(r).aggregateSingleCellFeatures(rec_array(r).Units);
    rec_array(r).saveObject();
end

%%
rec_array = recording_array_from_single_files(sorting_path_list(1:50));
for r = 1:length(rec_array)
    if isempty(rec_array(r).Units)
        disp(r)
            rec_array(r).Parameters.Save.Overwrite = true;
            rec_array(r).performAnalyses();
            rec_array(r).saveObject();
    end
end