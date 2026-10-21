function output_path_list = split_sortings(sorting_path, split_times, output_folder)
% SPLIT_SORTINGS  Split a spike sorting output into time-windowed sub-sortings.
%
% Reads spike_times.npy and spike_templates.npy (or spike_clusters.npy) from sorting_path,
% slices them into the requested windows, and writes each window to a new subfolder.
% Skips processing if the output files already exist.
%
% Each window gets its own params.py (n_samples = the window length, so SpikeData.Duration
% and firing rates refer to the window) and a split_window.txt marker, by which SpikeData
% gives all windows of one sorting the SortingID of that sorting (unit IDs are then shared
% between windows).
%
% INPUTS:
%   sorting_path  - Path to the root sorting folder (containing spike_times.npy etc.)
%   split_times   - N×2 double matrix of [start, end] times in minutes per window
%   output_folder - (optional) Base name for output subfolders; default "split"
%
% OUTPUTS:
%   output_path_list - String array of paths to the created sub-sorting folders

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

if ~all(arrayfun(@(x) exist(x,'file'),fullfile(output_path_list, "spike_times.npy")))

    spike_times = readNPY(fullfile(sorting_path, "spike_times.npy"));
    label_files = ["spike_templates.npy", "spike_clusters.npy", "amplitudes.npy"];
    label_files = label_files(arrayfun(@(f) isfile(fullfile(sorting_path, f)), label_files));
    assert(any(ismember(label_files, ["spike_templates.npy", "spike_clusters.npy"])), ...
        'split_sortings:noLabels', 'No spike_templates.npy or spike_clusters.npy in %s.', sorting_path);
    labels = arrayfun(@(f) readNPY(fullfile(sorting_path, f)), label_files, 'UniformOutput', false);

    params_txt = fileread(fullfile(sorting_path, "params.py"));
    tok = regexp(params_txt, '(?m)^\s*sample_rate\s*=\s*([-+0-9.eE]+)', 'tokens', 'once');
    assert(~isempty(tok), 'split_sortings:noSampleRate', 'No sample_rate in %s.', ...
        fullfile(sorting_path, "params.py"));
    sampling_rate = str2double(tok{1});

    window_samples = split_times * 60 * sampling_rate;

    for i = 1:size(window_samples,1)
        output_path = output_path_list(i);
        if ~isfolder(output_path), mkdir(output_path); end

        copyfile(fullfile(sorting_path, "channel_positions.npy"), output_path)
        copyfile(fullfile(sorting_path, "templates.npy"), output_path)

        % params.py with the window length as n_samples (the parent value would make the
        % duration, and every firing rate, refer to the whole parent recording)
        n_window = round(window_samples(i,2) - window_samples(i,1));
        txt = regexprep(params_txt, '(?m)^(\s*n_samples\w*\s*=\s*)[-+0-9.eE]+', sprintf('$1%d', n_window));
        if isempty(regexp(txt, '(?m)^\s*n_samples\w*\s*=', 'once'))
            txt = sprintf('%s%sn_samples = %d\n', txt, newline, n_window);
        end
        fid = fopen(fullfile(output_path, "params.py"), 'w');
        assert(fid > 0, 'split_sortings:write', 'Cannot write params.py in %s.', output_path);
        fwrite(fid, txt);
        fclose(fid);

        in_window = spike_times >= window_samples(i,1) & spike_times < window_samples(i,2);
        writeNPY(spike_times(in_window) - cast(window_samples(i,1), 'like', spike_times), ...
            fullfile(output_path, "spike_times.npy"))
        for k = 1:numel(label_files)
            writeNPY(labels{k}(in_window), fullfile(output_path, label_files(k)))
        end

        fid = fopen(fullfile(output_path, "split_window.txt"), 'w');
        fprintf(fid, 'parent = %s\nstart_min = %g\nend_min = %g\n', sorting_path, split_times(i,1), split_times(i,2));
        fclose(fid);
    end
else
    disp("Skipped splitting, output files already exist")
end
end
