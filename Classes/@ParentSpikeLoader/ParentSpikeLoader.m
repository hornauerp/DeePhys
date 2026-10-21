classdef ParentSpikeLoader
% PARENTSPIKELOADER  Shared session cache for parent-recording ACG computation.
%
% Static utility — no instantiation. For each unique combination of parent
% Kilosort directory + ACG parameters, computes a batch CCG once for ALL
% templates in that parent recording, extracts the auto-correlogram diagonal,
% and caches the result in memory.
%
% When multiple child RecordingProcessors share the same parent, this ensures
% the CCG is computed only once and each child simply looks up its own
% TemplateIDs from the returned map.
%
% USAGE:
%   params = struct('BinSize', 0.0005, 'Lag', 0.1);
%   acg_map = ParentSpikeLoader.loadFullACGs('/path/to/parent_ks', 30000, params);
%   acg = acg_map('42');   % ACG for TemplateID 42 (key is string)
%
%   ParentSpikeLoader.clearCache();  % free memory between experiments

    methods (Static)

        function acg_map = loadFullACGs(parent_path, sampling_rate, acg_params)
        % LOADFULLACGS  Return (or compute and cache) full-recording ACGs.
        %
        % INPUTS:
        %   parent_path   - absolute path to parent Kilosort directory
        %   sampling_rate - sampling rate in Hz (from child SpikeData)
        %   acg_params    - struct with fields BinSize (s) and Lag (s)
        %
        % OUTPUT:
        %   acg_map - containers.Map from TemplateID string → (n_bins×1) double ACG
            arguments
                parent_path   (1,1) string
                sampling_rate (1,1) double
                acg_params    (1,1) struct
            end

            cache = ParentSpikeLoader.getCache();
            key   = paramHash(struct( ...
                'ParentPath', char(parent_path), ...
                'BinSize',    acg_params.BinSize, ...
                'Lag',        acg_params.Lag));

            if cache.isKey(key)
                acg_map = cache(key);
                return
            end

            % Load full parent spike train
            npy_times = fullfile(parent_path, 'spike_times.npy');
            npy_units = fullfile(parent_path, 'spike_templates.npy');
            if ~isfile(npy_times) || ~isfile(npy_units)
                error('ParentSpikeLoader:missingFiles', ...
                    'spike_times.npy or spike_templates.npy not found in: %s', parent_path);
            end

            raw_times  = double(readNPY(npy_times));
            raw_units  = double(readNPY(npy_units));
            spike_times = raw_times / sampling_rate;
            % Map to 1-based TemplateIDs (SpikeData uses Kilosort value + 1). Kilosort writes
            % 0-based values, but parent exports such as qc_output are already 1-based; adding
            % 1 there gave every unit the ACG of the template before it.
            base = ParentSpikeLoader.templateIndexBase(parent_path, raw_units);
            spike_units = raw_units + (1 - base);

            % Sort by time (required by CCG)
            [spike_times, sort_idx] = sort(spike_times);
            spike_units = spike_units(sort_idx);

            % Compact local IDs (findgroups preserves ascending template order)
            [unit_ids_local, group_tids] = findgroups(spike_units);

            fprintf('ParentSpikeLoader: computing ACGs for %d templates in %s ...\n', ...
                numel(group_tids), parent_path);

            % One autocorrelogram per template from its own spike train. (The full
            % template x template CCG held n_templates^2 x n_bins counts only to read
            % its diagonal: ~3 GB for 1,000 templates at 401 bins.)
            train_of = splitapply(@(t) {t}, spike_times, unit_ids_local);
            acg_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for t = 1:numel(group_tids)
                raw_acg = CellTypeClassifier.computeACG(train_of{t}, acg_params.BinSize, ...
                    acg_params.Lag, sampling_rate);
                mx = max(raw_acg);
                if mx > 0
                    raw_acg = raw_acg / mx;
                end
                acg_map(num2str(group_tids(t))) = raw_acg;
            end

            cache(key) = acg_map;  %#ok<NASGU>
            fprintf('ParentSpikeLoader: done (%d templates cached).\n', numel(group_tids));
        end

        function clearCache()
        % CLEARCACHE  Remove all cached ACG maps to free memory.
            cache = ParentSpikeLoader.getCache();
            remove(cache, keys(cache));
        end

        function base = templateIndexBase(parent_path, raw_units)
        % TEMPLATEINDEXBASE  0 if spike_templates.npy is 0-based (Kilosort), 1 if 1-based.
        %
        % A 0 value means 0-based. Otherwise the values are compared with the number of
        % templates (templates.npy header, or phy_ids.npy for QC exports): a maximum equal to
        % that count is only possible with 1-based values. If neither file settles it, the
        % Kilosort convention is assumed with a warning.
            if any(raw_units == 0)
                base = 0;
                return
            end
            n_templates = NaN;
            tf = fullfile(parent_path, 'templates.npy');
            pf = fullfile(parent_path, 'phy_ids.npy');
            if isfile(tf)
                shape = readNPYheader(tf);
                n_templates = shape(1);
            elseif isfile(pf)
                n_templates = numel(readNPY(pf));
            end
            if ~isnan(n_templates) && max(raw_units) >= n_templates
                base = 1;
            else
                base = 0;
                if isnan(n_templates)
                    warning('ParentSpikeLoader:indexBaseAmbiguous', ...
                        ['%s: spike_templates has no 0 and no templates.npy/phy_ids.npy to check ' ...
                         'against; assuming 0-based (Kilosort) template values.'], parent_path);
                else
                    % Values 1..n-1 of n templates: 0-based with template 0 silent, or 1-based
                    % with the last template silent. Both fit; a wrong guess shifts every unit's
                    % parent ACG to the neighbouring template.
                    warning('ParentSpikeLoader:indexBaseAmbiguous', ...
                        ['%s: spike_templates has no 0 and its maximum (%d) is below the template ' ...
                         'count (%d), so 0- and 1-based numbering both fit; assuming 0-based ' ...
                         '(Kilosort). Check the parent ACGs if this sorting was exported 1-based.'], ...
                        parent_path, max(raw_units), n_templates);
                end
            end
        end

    end

    methods (Static, Access = private)

        function cache = getCache()
        % Persistent session-scoped cache (containers.Map of ACG maps).
            persistent c
            if isempty(c)
                c = containers.Map('KeyType', 'char', 'ValueType', 'any');
            end
            cache = c;
        end

    end
end
