function sorting_path = make_mock_sorting(output_dir, n_units, n_spikes_per_unit, sampling_rate, duration_s)
% MAKE_MOCK_SORTING  Create a minimal Kilosort-compatible output directory for testing.
%
% Generates the minimum set of files that MEArecording needs to load a
% recording (spike_times.npy, spike_templates.npy, templates.npy,
% channel_positions.npy, params.py). All waveforms are random noise with a
% plausible negative deflection to pass basic waveform QC.
%
% INPUTS:
%   output_dir          - Parent directory in which to create the mock folder
%   n_units             - Number of units to simulate          (default 5)
%   n_spikes_per_unit   - Spikes per unit                      (default 50)
%   sampling_rate       - Sampling rate in Hz                  (default 20000)
%   duration_s          - Recording duration in seconds        (default 60)
%
% OUTPUT:
%   sorting_path - Full path to the created sorting directory
%
% EXAMPLE:
%   path = make_mock_sorting(tempdir);
%   metadata.InputPath = path;
%   obj = MEArecording(metadata);   % smoke-test full construction

arguments
    output_dir          (1,1) string
    n_units             (1,1) double = 5
    n_spikes_per_unit   (1,1) double = 50
    sampling_rate       (1,1) double = 20000
    duration_s          (1,1) double = 60
end

sorting_path = fullfile(output_dir, 'mock_sorting');
if ~isfolder(sorting_path)
    mkdir(sorting_path);
end

rng(0);   % reproducible waveform noise

% ------------------------------------------------------------------ %
%  spike_times.npy  (N_spikes × 1, uint64, sample indices, sorted)
%  spike_templates.npy  (N_spikes × 1, uint32, 0-indexed unit IDs)
% ------------------------------------------------------------------ %
n_spikes     = n_units * n_spikes_per_unit;
max_sample   = floor(sampling_rate * duration_s);
all_times    = sort(randi([sampling_rate, max_sample], n_spikes, 1), 'ascend');
templates_id = repelem((0 : n_units-1)', n_spikes_per_unit);   % 0-indexed

writeNPY(uint64(all_times),  fullfile(sorting_path, 'spike_times.npy'));
writeNPY(uint32(templates_id), fullfile(sorting_path, 'spike_templates.npy'));

% ------------------------------------------------------------------ %
%  templates.npy  (n_units × n_samples × n_channels, single)
%  82 samples ≈ 4.1 ms at 20 kHz; 4 channels
% ------------------------------------------------------------------ %
n_samples  = 82;
n_channels = 4;
waveforms  = zeros(n_units, n_samples, n_channels, 'single');

% Smooth biphasic spike (trough near sample 20, slower positive repolarisation), decaying
% across channels, plus a little noise. Amplitudes are scaled by the QC LSB (6.2 uV) and
% must clear RecordingProcessor's default Amplitude window [20 1000]; white-noise
% templates would fail the waveform noise check (too many zero crossings).
t = 1:n_samples;
for u = 1:n_units
    amp   = 5 + 3 * rand;
    shape = -amp * exp(-(t - 20).^2 / (2 * 1.5^2)) + 0.3 * amp * exp(-(t - 30).^2 / (2 * 4^2));
    for c = 1:n_channels
        waveforms(u, :, c) = single(shape * 0.6^(c - 1) + 0.002 * amp * randn(1, n_samples));
    end
end

writeNPY(waveforms, fullfile(sorting_path, 'templates.npy'));

% ------------------------------------------------------------------ %
%  channel_positions.npy  (n_channels × 2, single)
%  Simple linear probe with 17.32 µm pitch
% ------------------------------------------------------------------ %
positions = single([(0:n_channels-1)' * 17.32, zeros(n_channels, 1)]);
writeNPY(positions, fullfile(sorting_path, 'channel_positions.npy'));

% ------------------------------------------------------------------ %
%  params.py
% ------------------------------------------------------------------ %
fid = fopen(fullfile(sorting_path, 'params.py'), 'w');
fprintf(fid, 'sample_rate = %d\n',    sampling_rate);
fprintf(fid, 'n_channels_dat = %d\n', n_channels);
fprintf(fid, 'dtype = ''int16''\n');
fclose(fid);

fprintf('Mock sorting written to: %s\n  %i units, %i spikes/unit, %i Hz, %i s\n', ...
    sorting_path, n_units, n_spikes_per_unit, sampling_rate, duration_s);
end
