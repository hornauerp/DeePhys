function [is_valid, issues] = validate(obj)
% VALIDATE  Check MEArecording for common data integrity issues.
%
%   [is_valid, issues] = obj.validate()
%
%   Returns true if no critical issues found. Issues is a struct array
%   with fields: severity ("error"|"warning"), message, field.
%
%   Can be called at any point — in the constructor (after parseRecordingInfo),
%   after loading from file, or manually for QC.

    issues = struct('severity', {}, 'message', {}, 'field', {});

    % --- Critical checks (would cause downstream failures) ---

    % Recording duration
    if ~isfield(obj.RecordingInfo, 'Duration') || isempty(obj.RecordingInfo.Duration)
        issues(end+1) = struct('severity', 'error', 'message', 'RecordingInfo.Duration is empty.', 'field', 'RecordingInfo.Duration');
    elseif obj.RecordingInfo.Duration <= 0
        issues(end+1) = struct('severity', 'error', 'message', sprintf('RecordingInfo.Duration is non-positive: %.4f.', obj.RecordingInfo.Duration), 'field', 'RecordingInfo.Duration');
    end

    % Sampling rate
    if ~isfield(obj.RecordingInfo, 'SamplingRate') || isempty(obj.RecordingInfo.SamplingRate)
        issues(end+1) = struct('severity', 'error', 'message', 'RecordingInfo.SamplingRate is empty.', 'field', 'RecordingInfo.SamplingRate');
    elseif obj.RecordingInfo.SamplingRate <= 0
        issues(end+1) = struct('severity', 'error', 'message', sprintf('RecordingInfo.SamplingRate is non-positive: %.1f.', obj.RecordingInfo.SamplingRate), 'field', 'RecordingInfo.SamplingRate');
    end

    % Spike times
    if isfield(obj.Spikes, 'Times')
        st = obj.Spikes.Times;
        if isempty(st)
            issues(end+1) = struct('severity', 'warning', 'message', 'No spike times found.', 'field', 'Spikes.Times');
        elseif any(isnan(st))
            issues(end+1) = struct('severity', 'error', 'message', sprintf('%d NaN values in spike times.', sum(isnan(st))), 'field', 'Spikes.Times');
        elseif any(st < 0)
            issues(end+1) = struct('severity', 'error', 'message', sprintf('%d negative spike times found.', sum(st < 0)), 'field', 'Spikes.Times');
        end
    end

    % Spike templates
    if isfield(obj.Spikes, 'Templates')
        templates = obj.Spikes.Templates;
        if isempty(templates)
            issues(end+1) = struct('severity', 'warning', 'message', 'No spike templates found.', 'field', 'Spikes.Templates');
        end
    end

    % --- Unit checks ---
    if ~isempty(obj.Units)
        n_units = length(obj.Units);

        % Check for units with no spikes
        empty_units = 0;
        for u = 1:n_units
            if isempty(obj.Units(u).SpikeTimes)
                empty_units = empty_units + 1;
            end
        end
        if empty_units > 0
            issues(end+1) = struct('severity', 'warning', 'message', ...
                sprintf('%d / %d units have empty spike times.', empty_units, n_units), 'field', 'Units');
        end

        % Check for duplicate template IDs
        template_ids = [obj.Units.TemplateID];
        if length(unique(template_ids)) < length(template_ids)
            issues(end+1) = struct('severity', 'warning', 'message', ...
                'Duplicate TemplateIDs found among units.', 'field', 'Units.TemplateID');
        end

        % Check for NaN in waveforms
        nan_wf = 0;
        for u = 1:n_units
            if any(isnan(obj.Units(u).ReferenceWaveform), 'all')
                nan_wf = nan_wf + 1;
            end
        end
        if nan_wf > 0
            issues(end+1) = struct('severity', 'warning', 'message', ...
                sprintf('%d / %d units have NaN in ReferenceWaveform.', nan_wf, n_units), 'field', 'Units.ReferenceWaveform');
        end
    end

    % --- Metadata checks ---
    if ~isempty(obj.Metadata)
        if ~isfield(obj.Metadata, 'InputPath') || isempty(obj.Metadata.InputPath)
            issues(end+1) = struct('severity', 'error', 'message', 'Metadata.InputPath is empty.', 'field', 'Metadata.InputPath');
        end
    else
        issues(end+1) = struct('severity', 'error', 'message', 'Metadata is empty.', 'field', 'Metadata');
    end

    % --- Bombcell checks ---
    if ~isempty(obj.Units) && ~isempty(obj.Units(1).BombcellMetrics)
        n_units = length(obj.Units);
        bc_table = vertcat(obj.Units.BombcellMetrics);

        % Check for units where bombcell metrics are all NaN (failed computation)
        all_nan_bc = sum(isnan(bc_table.BC_nPeaks) & isnan(bc_table.BC_nTroughs) & isnan(bc_table.BC_fractionRPVs));
        if all_nan_bc > 0
            issues(end+1) = struct('severity', 'warning', 'message', ...
                sprintf('%d / %d units have all-NaN bombcell metrics (computation may have failed).', all_nan_bc, n_units), ...
                'field', 'Units.BombcellMetrics');
        end

        % Check for high RPV contamination among accepted units
        rpv_vals = bc_table.BC_fractionRPVs;
        high_rpv = sum(rpv_vals > 0.1 & ~isnan(rpv_vals));
        if high_rpv > 0
            issues(end+1) = struct('severity', 'warning', 'message', ...
                sprintf('%d / %d units have bombcell RPV > 10%% (potential contamination).', high_rpv, n_units), ...
                'field', 'Units.BombcellMetrics.BC_fractionRPVs');
        end

        % Check for low presence ratio
        pr_vals = bc_table.BC_presenceRatio;
        low_pr = sum(pr_vals < 0.5 & ~isnan(pr_vals));
        if low_pr > 0
            issues(end+1) = struct('severity', 'warning', 'message', ...
                sprintf('%d / %d units have bombcell presence ratio < 0.5 (intermittent units).', low_pr, n_units), ...
                'field', 'Units.BombcellMetrics.BC_presenceRatio');
        end

        % Check for mismatched BombcellType
        bc_types = [obj.Units.BombcellType];
        if any(isnan(bc_types))
            issues(end+1) = struct('severity', 'warning', 'message', ...
                sprintf('%d / %d units have NaN BombcellType (classification failed).', sum(isnan(bc_types)), n_units), ...
                'field', 'Units.BombcellType');
        end
    end

    % --- Connectivity checks ---
    if ~isempty(obj.Connectivity) && isfield(obj.Connectivity, 'CCG')
        ccg = obj.Connectivity.CCG;
        if isfield(ccg, 'CCGs') && any(isnan(ccg.CCGs), 'all')
            issues(end+1) = struct('severity', 'warning', 'message', 'NaN values in CCG matrix.', 'field', 'Connectivity.CCG');
        end
    end

    % Determine overall validity
    severities = {issues.severity};
    is_valid = ~any(strcmp(severities, 'error'));

    % Print summary
    n_errors = sum(strcmp(severities, 'error'));
    n_warnings = sum(strcmp(severities, 'warning'));
    if n_errors > 0
        fprintf('Validation FAILED: %d error(s), %d warning(s)\n', n_errors, n_warnings);
        for i = 1:length(issues)
            if strcmp(issues(i).severity, 'error')
                fprintf('  ERROR: %s\n', issues(i).message);
            end
        end
    elseif n_warnings > 0
        fprintf('Validation passed with %d warning(s)\n', n_warnings);
    end
end
