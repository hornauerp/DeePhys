       function unit_array = generateUnits(obj)
           [max_amplitudes, reference_electrode, norm_wf_matrix] = obj.generateWaveformMatrix();
           [firing_rates, unit_spike_times] = obj.calculateFiringRates(length(max_amplitudes));

           % Run bombcell QC on all templates (before filtering)
           run_bombcell = isfield(obj.Parameters.QC, 'Bombcell') && obj.Parameters.QC.Bombcell.Enable;
           if run_bombcell
               [bombcell_metrics, bombcell_types] = obj.runBombcellQC(reference_electrode);
           end

           % Check if KS label filtering is requested
           use_ks = isfield(obj.Parameters.QC, 'KSLabel') && obj.Parameters.QC.KSLabel.Enable;
           ks_labels = string.empty; ks_ids = [];
           if use_ks
               try
                   [ks_labels, ks_ids] = obj.getKSLabels();
               catch ME
                   warning("DeePhys:ksLabelLoad", "Could not load KSLabels: %s", ME.message);
                   use_ks = false;
               end
           end

           if isempty(obj.Parameters.QC.GoodUnits)
               good_units = performUnitQC(obj, max_amplitudes, firing_rates', unit_spike_times', norm_wf_matrix);

               % Apply Kilosort label filter
               if use_ks
                   n_templates = length(max_amplitudes);
                   ks_pass = false(n_templates, 1);
                   accepted_labels = string(obj.Parameters.QC.KSLabel.AcceptedLabels);
                   for k = 1:length(ks_ids)
                       tid = ks_ids(k) + 1; % 0-indexed -> 1-indexed
                       if tid >= 1 && tid <= n_templates
                           ks_pass(tid) = ismember(ks_labels(k), accepted_labels);
                       end
                   end
                   n_before = sum(good_units);
                   good_units = good_units & ks_pass;
                   fprintf('KSLabel filter kept %s units: %d -> %d\n', ...
                       strjoin(accepted_labels, "/"), n_before, sum(good_units))
               end

               % Apply bombcell filter: reject units bombcell classifies as unacceptable
               if run_bombcell && obj.Parameters.QC.Bombcell.FilterUnits
                   accepted = obj.Parameters.QC.Bombcell.AcceptedTypes;
                   bombcell_pass = ismember(bombcell_types, accepted);
                   n_before = sum(good_units);
                   good_units = good_units & bombcell_pass(:);
                   fprintf('Bombcell filtered %d additional units (%d -> %d)\n', ...
                       n_before - sum(good_units), n_before, sum(good_units))
               end

               obj.Parameters.QC.GoodUnits = good_units;
           else
               good_units = obj.Parameters.QC.GoodUnits;
               if length(good_units) ~= length(firing_rates)
                   error("Good unit IDs do not match. Did you use the correct ones?")
               end
               fprintf('Kept the %i good units provided\n', sum(good_units))
           end

           if sum(good_units) >= obj.Parameters.QC.N_Units
               good_amplitudes = max_amplitudes(good_units);
               good_unit_spike_times = unit_spike_times(good_units);
               good_wf_matrix = norm_wf_matrix(:,good_units);
               waveform_features = obj.inferWaveformFeatures(good_amplitudes, good_wf_matrix);
               good_template_ids = find(good_units);
               good_reference_electrodes = reference_electrode(good_units);
               unit_array = Unit();
               for u = 1:length(waveform_features)
                   unit_array(u) = Unit(obj, good_wf_matrix(:,u), good_reference_electrodes(u), good_unit_spike_times{u}, waveform_features{u});
                   unit_array(u).TemplateID = good_template_ids(u);
               end

               % Store Kilosort labels on each unit
               if use_ks
                   for u = 1:length(unit_array)
                       tid_0 = good_template_ids(u) - 1; % 1-indexed -> 0-indexed for lookup
                       match = find(ks_ids == tid_0, 1);
                       if ~isempty(match)
                           unit_array(u).KSLabel = ks_labels(match);
                       end
                   end
               end

               % Store bombcell metrics as feature table on each unit
               if run_bombcell
                   bc_names = Unit.returnFeatureNames("bc");
                   for u = 1:length(unit_array)
                       tid = good_template_ids(u);
                       bc_vals = [bombcell_metrics.nPeaks(tid), ...
                                  bombcell_metrics.nTroughs(tid), ...
                                  bombcell_metrics.waveformDuration_peakTrough(tid), ...
                                  bombcell_metrics.mainPeakToTroughRatio(tid), ...
                                  bombcell_metrics.peak1ToPeak2Ratio(tid), ...
                                  bombcell_metrics.scndPeakToTroughRatio(tid), ...
                                  bombcell_metrics.troughToPeak2Ratio(tid), ...
                                  bombcell_metrics.waveformBaselineFlatness(tid), ...
                                  bombcell_metrics.spatialDecaySlope(tid), ...
                                  bombcell_metrics.fractionRPVs_estimatedTauR(tid), ...
                                  bombcell_metrics.presenceRatio(tid), ...
                                  bombcell_metrics.percentageSpikesMissing_gaussian(tid)];
                       unit_array(u).BombcellMetrics = array2table(bc_vals, 'VariableNames', bc_names);
                       unit_array(u).BombcellType = bombcell_types(tid);
                   end
               end

               no_act_idx = arrayfun(@(x) isempty(x.ActivityFeatures),unit_array);
               %If activity features were not computed (too few spikes), fill with zeros to prevent issues when tracking units
               act_features = Unit.returnFeatureNames("act");
                empty_act_table = array2table(zeros(1,length(act_features)),'VariableNames',act_features);
               [unit_array(no_act_idx).ActivityFeatures] = deal(empty_act_table);
               if obj.Parameters.Analyses.Regularity
                   no_reg_idx = arrayfun(@(x) isempty(x.RegularityFeatures),unit_array);
                   reg_features = Unit.returnFeatureNames("reg");
                   empty_reg_table = array2table(zeros(1,length(reg_features)),'VariableNames',reg_features);
                   [unit_array(no_reg_idx).RegularityFeatures] = deal(empty_reg_table);
               end
               if obj.Parameters.Analyses.Catch22
                   no_c22_idx = arrayfun(@(x) isempty(x.Catch22),unit_array);
                   c22_features = Unit.returnFeatureNames("c22");
                   empty_c22_table = array2table(zeros(1,length(c22_features)),'VariableNames',"SC_" + c22_features);
                   [unit_array(no_c22_idx).Catch22] = deal(empty_c22_table);
               end

               obj.Units = unit_array;
               obj.updateSpikeTimes();
           else
               warning("Not enough good units found (%d passed QC, %d required)", ...
                   sum(good_units), obj.Parameters.QC.N_Units)
               unit_array = Unit.empty(0,0);
           end
       end
