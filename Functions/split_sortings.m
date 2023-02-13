function output_path_list = split_sortings(sorting_path, split_times, output_folder)
arguments
    sorting_path string
    split_times (:,2) double %In minutes; first column indicates start times, second column corresponding end times
    output_folder string = "split"
end

if length(output_folder) < size(split_times,1)
    output_folders = output_folder + "_" + split_times(:,1) + "-" + split_times(:,2);
else
    output_folders = output_folder;
end

output_path_list = fullfile(sorting_path, output_folders);

if ~all(arrayfun(@(x) exist(x,'file'),fullfile(output_path_list, "spike_templates.npy")))
    
    spike_times_file = fullfile(sorting_path, "spike_times.npy");
    spike_templates_file = fullfile(sorting_path, "spike_templates.npy");
    
    spike_times = readNPY(spike_times_file);
    spike_templates = readNPY(spike_templates_file);
    
    params_path = fullfile(sorting_path, "params.py");
    fileID = fopen(params_path,'r');
    params_txt = fscanf(fileID,'%c');
    params_parts = strsplit(params_txt,'\n');
    sr_idx = startsWith(params_parts,'sample');
    sr = strsplit(params_parts{sr_idx},' = ');
    sampling_rate = str2double(sr{2});
    
    split_times = split_times * 60 * sampling_rate;
    
    for i = 1:size(split_times,1)
        output_path = output_path_list(i);
        mkdir(output_path)
        
        copyfile(fullfile(sorting_path, "channel_positions.npy"), output_path)
        copyfile(fullfile(sorting_path, "templates.npy"), output_path)
        copyfile(fullfile(sorting_path, "params.py"), output_path)
        
        split_spike_times = spike_times(spike_times > split_times(i,1) & spike_times < split_times(i,2)) - split_times(i,1);
        split_spike_templates = spike_templates(spike_times > split_times(i,1) & spike_times < split_times(i,2));
        
        writeNPY(split_spike_times, fullfile(output_path, "spike_times.npy"))
        writeNPY(split_spike_templates, fullfile(output_path, "spike_templates.npy"))
        
    end
else
    disp("Output files already exist")
end
end

