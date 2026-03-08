       function unit_array = generateUnits(obj)
           [max_amplitudes, reference_electrode, norm_wf_matrix] = obj.generateWaveformMatrix();
           [firing_rates, unit_spike_times] = obj.calculateFiringRates(length(max_amplitudes));


           if isempty(obj.Parameters.QC.GoodUnits)
               good_units = performUnitQC(obj, max_amplitudes, firing_rates', unit_spike_times', norm_wf_matrix);
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
               warning("Not enough good units found")
           end
       end
