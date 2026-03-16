function [qMetric, unitType] = runBombcellQC(obj, reference_electrodes)
% RUNBOMBCELLQC  Compute bombcell quality metrics using spike-sorted templates.
%
%   [qMetric, unitType] = obj.runBombcellQC(reference_electrodes)
%
%   Runs bombcell quality metrics on all templates without requiring raw
%   recordings. Uses template waveforms for morphology analysis and spike
%   times for temporal quality checks.
%
%   Computed metrics:
%     - Refractory period violations (Hill et al. method)
%     - Presence ratio
%     - Percentage spikes missing (requires amplitudes.npy)
%     - Spike count
%     - Waveform shape (nPeaks, nTroughs, duration, ratios, baseline)
%     - Spatial decay slope (if channel geometry permits)
%
%   Skipped (require raw recording or PC features):
%     - rawAmplitude, signalToNoiseRatio
%     - isolationDistance, Lratio, silhouetteScore
%     - maxDriftEstimate, cumDriftEstimate
%
%   Inputs:
%     reference_electrodes - [nTemplates x 1] peak channel index per template
%
%   Outputs:
%     qMetric  - struct with [nTemplates x 1] fields for each metric
%     unitType - [nTemplates x 1] bombcell classification:
%                0 = noise, 1 = good, 2 = MUA,
%                3 = non-somatic good, 4 = non-somatic MUA

    bc_cfg = obj.Parameters.QC.Bombcell;

    % Add bombcell to path if specified
    if isfield(bc_cfg, 'Path') && ~isempty(bc_cfg.Path)
        addpath(genpath(bc_cfg.Path));
    end

    % Load original template waveforms (before DeePhys zero-column removal)
    templateWaveforms = obj.getTemplateMatrix();
    n_units = size(templateWaveforms, 1);
    spike_width = size(templateWaveforms, 2);

    % Load spike data
    [spike_times, spike_units] = obj.getSpikeTimes();

    % Load per-spike amplitudes if available
    amp_path = fullfile(obj.Metadata.InputPath, 'amplitudes.npy');
    has_amplitudes = exist(amp_path, 'file') == 2;
    if has_amplitudes
        all_amplitudes = double(readNPY(amp_path));
        disp("Loaded per-spike amplitudes for bombcell QC")
    else
        warning("DeePhys:noBombcellAmplitudes", ...
            "amplitudes.npy not found — percentageSpikesMissing will be unavailable")
    end

    % Build bombcell parameter struct
    param = buildBombcellParam(obj, spike_width);

    % Channel positions for spatial decay
    channelPositions = obj.RecordingInfo.ElectrodeCoordinates;

    % Baseline window (sample indices)
    waveformBaselineWindow = [param.waveformBaselineWindowStart, ...
                              param.waveformBaselineWindowStop];

    duration = obj.RecordingInfo.Duration;
    timeChunks = [0, duration];

    % --- Pre-allocate metric fields ---
    metric_fields = {'nSpikes', 'fractionRPVs_estimatedTauR', 'presenceRatio', ...
        'percentageSpikesMissing_gaussian', 'percentageSpikesMissing_symmetric', ...
        'ksTest_pValue', 'nPeaks', 'nTroughs', 'spatialDecaySlope', ...
        'waveformBaselineFlatness', 'waveformDuration_peakTrough', ...
        'mainPeakToTroughRatio', 'peak1ToPeak2Ratio', 'scndPeakToTroughRatio', ...
        'troughToPeak2Ratio', 'mainPeak_before_width', 'mainPeak_after_width', ...
        'mainTrough_width'};
    for f = 1:length(metric_fields)
        qMetric.(metric_fields{f}) = NaN(n_units, 1);
    end

    % --- Compute per-unit metrics ---
    for u = 1:n_units
        unit_mask = spike_units == u;
        theseSpikeTimes = spike_times(unit_mask);
        n_spikes = length(theseSpikeTimes);

        % Spike count (always available)
        qMetric.nSpikes(u) = n_spikes;

        if n_spikes < 2
            continue
        end

        % Per-spike amplitudes
        if has_amplitudes
            theseAmplitudes = all_amplitudes(unit_mask);
        else
            theseAmplitudes = ones(n_spikes, 1);
        end

        % --- Refractory period violations (Hill method) ---
        tauR = param.tauR_valuesMin;
        try
            [RPVrate, ~, ~] = bc.qm.fractionRPviolations(...
                theseSpikeTimes, theseAmplitudes, tauR, param, timeChunks);
            qMetric.fractionRPVs_estimatedTauR(u) = RPVrate(1);
        catch ME
            warning("DeePhys:bombcellRPV", "RPV failed for unit %d: %s", u, ME.message);
        end

        % --- Presence ratio ---
        try
            qMetric.presenceRatio(u) = bc.qm.presenceRatio(...
                theseSpikeTimes, theseAmplitudes, ...
                param.presenceRatioBinSize, 0, duration, param);
        catch ME
            warning("DeePhys:bombcellPresence", "Presence ratio failed for unit %d: %s", u, ME.message);
        end

        % --- Percentage spikes missing (needs real amplitudes) ---
        if has_amplitudes && n_spikes >= 20
            try
                [pMissGauss, pMissSym, ksPval, ~, ~, ~] = bc.qm.percSpikesMissing(...
                    theseAmplitudes, theseSpikeTimes, timeChunks, param);
                qMetric.percentageSpikesMissing_gaussian(u) = pMissGauss(1);
                qMetric.percentageSpikesMissing_symmetric(u) = pMissSym(1);
                qMetric.ksTest_pValue(u) = ksPval(1);
            catch ME
                warning("DeePhys:bombcellMissing", "percSpikesMissing failed for unit %d: %s", u, ME.message);
            end
        end

        % --- Waveform shape ---
        try
            [nPeaks, nTroughs, spatialDecaySlope, waveformBaseline, ...
             scndPeakToTroughRatio, mainPeakToTroughRatio, peak1ToPeak2Ratio, ...
             troughToPeak2Ratio, mainPeak_before_width, mainPeak_after_width, ...
             mainTrough_width, ~, ~, waveformDuration_peakTrough] = ...
                bc.qm.waveformShape(templateWaveforms, u, reference_electrodes(u), ...
                param, channelPositions, waveformBaselineWindow);

            qMetric.nPeaks(u) = nPeaks;
            qMetric.nTroughs(u) = nTroughs;
            qMetric.spatialDecaySlope(u) = spatialDecaySlope;
            qMetric.waveformBaselineFlatness(u) = waveformBaseline;
            qMetric.waveformDuration_peakTrough(u) = waveformDuration_peakTrough;
            qMetric.mainPeakToTroughRatio(u) = mainPeakToTroughRatio;
            qMetric.peak1ToPeak2Ratio(u) = peak1ToPeak2Ratio;
            qMetric.scndPeakToTroughRatio(u) = scndPeakToTroughRatio;
            qMetric.troughToPeak2Ratio(u) = troughToPeak2Ratio;
            qMetric.mainPeak_before_width(u) = mainPeak_before_width;
            qMetric.mainPeak_after_width(u) = mainPeak_after_width;
            qMetric.mainTrough_width(u) = mainTrough_width;
        catch ME
            warning("DeePhys:bombcellWaveform", "Waveform shape failed for unit %d: %s", u, ME.message);
        end
    end

    % --- Add fields expected by getQualityUnitType ---
    qMetric.phy_clusterID = (0:n_units-1)';
    qMetric.clusterID = (1:n_units)';
    qMetric.maxChannels = reference_electrodes(:);

    % Placeholders for unavailable metrics (raw recording / PC features)
    qMetric.rawAmplitude = NaN(n_units, 1);
    qMetric.signalToNoiseRatio = NaN(n_units, 1);
    qMetric.isolationDistance = NaN(n_units, 1);
    qMetric.Lratio = NaN(n_units, 1);
    qMetric.silhouetteScore = NaN(n_units, 1);
    qMetric.maxDriftEstimate = NaN(n_units, 1);
    qMetric.cumDriftEstimate = NaN(n_units, 1);
    qMetric.useTheseTimesStart = zeros(n_units, 1);
    qMetric.useTheseTimesStop = ones(n_units, 1) * duration;

    % --- Classify units ---
    try
        [unitType, ~] = bc.qm.getQualityUnitType(param, qMetric, tempdir());
    catch ME
        warning("DeePhys:bombcellClassify", ...
            "Unit classification failed: %s\nFalling back to NaN.", ME.message);
        unitType = NaN(n_units, 1);
    end

    fprintf('Bombcell QC: %d noise, %d good, %d MUA out of %d templates\n', ...
        sum(unitType == 0), sum(unitType == 1), sum(unitType == 2), n_units);
end


function param = buildBombcellParam(obj, spike_width)
% BUILDBOMBCELLPARAM  Construct bombcell-compatible parameter struct.
%   Builds a complete parameter struct for bombcell QC functions using
%   DeePhys recording info and bombcell default thresholds. User overrides
%   from Parameters.QC.Bombcell.Overrides are applied last.

    sr = obj.RecordingInfo.SamplingRate;
    bc_cfg = obj.Parameters.QC.Bombcell;

    % --- Core settings ---
    param.ephys_sample_rate = sr;
    param.spikeWidth = spike_width;
    param.plotDetails = 0;
    param.plotGlobal = 0;
    param.verbose = 0;

    % --- Disable features requiring raw data or PC features ---
    param.extractRaw = 0;
    param.computeDistanceMetrics = 0;
    param.computeTimeChunks = 0;
    param.computeDrift = 0;

    % --- Refractory period ---
    param.tauR_valuesMin = 0.002;   % 2 ms
    param.tauR_valuesStep = 0.0005;
    param.tauR_valuesMax = 0.002;   % Single value
    param.tauC = 0.0001;            % 0.1 ms censored period
    param.hillOrLlobetMethod = 1;   % Hill et al.

    % --- Presence ratio ---
    param.presenceRatioBinSize = 60; % seconds

    % --- Waveform analysis ---
    param.minThreshDetectPeaksTroughs = 0.2;
    param.normalizeSpDecay = 1;
    param.spDecayLinFit = 0;          % Exponential fit
    param.computeSpatialDecay = 1;

    % Baseline window scaled to actual spike width (82 is KS2/3 default)
    scale = spike_width / 82;
    param.waveformBaselineWindowStart = round(20 * scale);
    param.waveformBaselineWindowStop = round(30 * scale);
    param.maxWvBaselineFraction = 0.3;

    % --- Classification thresholds (noise) ---
    param.maxNPeaks = 2;
    param.maxNTroughs = 1;
    param.minWvDuration = 100;          % microseconds
    param.maxWvDuration = 1150;         % microseconds
    param.maxScndPeakToTroughRatio_noise = 0.8;
    param.minSpatialDecaySlope = -0.008;         % linear: a.u./um
    param.minSpatialDecaySlopeExp = 0.01;        % exponential
    param.maxSpatialDecaySlopeExp = 0.1;

    % --- Classification thresholds (MUA / good) ---
    param.maxPercSpikesMissing = 20;    % percent
    param.maxRPVviolations = 0.1;       % fraction
    param.minNumSpikes = 300;
    param.minPresenceRatio = 0.7;
    param.minAmplitude = 40;            % uV (unused: extractRaw=0)
    param.minSNR = 5;                   % unused: extractRaw=0

    % --- Non-somatic classification ---
    param.maxPeak1ToPeak2Ratio_nonSomatic = 3;
    param.maxMainPeakToTroughRatio_nonSomatic = 0.8;
    param.minWidthFirstPeak_nonSomatic = round(4 * scale);
    param.minWidthMainTrough_nonSomatic = round(5 * scale);
    param.minTroughToPeak2Ratio_nonSomatic = 5;
    param.splitGoodAndMua_NonSomatic = 0;

    % --- Apply user overrides ---
    if isfield(bc_cfg, 'Overrides') && isstruct(bc_cfg.Overrides) && ~isempty(fieldnames(bc_cfg.Overrides))
        fnames = fieldnames(bc_cfg.Overrides);
        for i = 1:length(fnames)
            param.(fnames{i}) = bc_cfg.Overrides.(fnames{i});
        end
    end
end
