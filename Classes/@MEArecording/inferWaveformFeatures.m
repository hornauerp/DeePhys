        function waveform_features = inferWaveformFeatures(obj, max_amplitudes, norm_wf_matrix)
            interpolation_factor = 10;
            ms_conversion = 10 * interpolation_factor; % 10 samples/ms * interpolation_factor = samples per ms in the interpolated waveform
            x = 1:size(norm_wf_matrix,1);
            xq = 1:(1/interpolation_factor):size(norm_wf_matrix,1);
            interp_wf_matrix = double(interp1(x,norm_wf_matrix,xq,'pchip'));
            zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %Function to detect zero crossings
            tx = zci(interp_wf_matrix); %Apply to interpolated waveform matrix
            [sample_idx,electrode_idx] = ind2sub(size(interp_wf_matrix),tx);
            unit_zero_crossings = splitapply(@(x) {x},sample_idx,electrode_idx);
            [unit_trough_value,unit_trough_idx] = min(interp_wf_matrix);
            peak_1_cutout = interp_wf_matrix(unit_trough_idx - ms_conversion:unit_trough_idx,:); %Limit detection range for peaks
            peak_2_cutout = interp_wf_matrix(unit_trough_idx:unit_trough_idx + ms_conversion,:);
            [peak_1_value, ~] = max(peak_1_cutout);
            [peak_2_value, peak_2_idx] = max(peak_2_cutout);
            peak_2_idx = peak_2_idx + unit_trough_idx - 1;
            asymmetry = (peak_2_value - peak_1_value) ./ (peak_2_value + peak_1_value);
            t2pdelay = (peak_2_idx - unit_trough_idx) ./ (obj.RecordingInfo.SamplingRate / 1000); %in [ms]
            t2pratio = abs(unit_trough_value ./ peak_2_value);
            waveform_features = cell(1,length(max_amplitudes));
            for u = 1:size(interp_wf_matrix,2)
                unit_features = struct();
                zc = [1 unit_zero_crossings{u}' length(xq)];
                [~,zc_pre_trough] = max(zc(zc<unit_trough_idx(u))); %Find zero crossing towards the trough
                zc_pre_trough = max([zc_pre_trough 2]);
                unit_features.AUC_peak_1 = trapz(interp_wf_matrix(zc(zc_pre_trough - 1):zc(zc_pre_trough)));
                unit_features.AUC_trough = abs(trapz(interp_wf_matrix(zc(zc_pre_trough):zc(zc_pre_trough + 1))));
                if length(zc(zc_pre_trough:end)) > 3 %Ensure that there is another zero crossing at the end
                    unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):zc(zc_pre_trough + 2)));
                else %Otherwise we integrate until the end of the waveform
                    unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):size(interp_wf_matrix,1)));
                end
                %We pad the signal with min/max values to ensure reliable slewrate
                %calculations
                rise_cutout = interp_wf_matrix(unit_trough_idx(u):peak_2_idx(u),u);
                padded_rise = [ones(100,1)*rise_cutout(1); rise_cutout; ones(100,1)*rise_cutout(end)];
                unit_features.Rise = mean(slewrate(padded_rise,interpolation_factor));
                decay_cutout = interp_wf_matrix(peak_2_idx(u):end,u);
                decay_cutout = decay_cutout(1:find(decay_cutout==min(decay_cutout), 1, 'first'));
                padded_decay = [ones(100,1)*decay_cutout(1); decay_cutout; ones(100,1)*decay_cutout(end)];
                unit_features.Decay = mean(slewrate(padded_decay,interpolation_factor));
                unit_features.Asymmetry = asymmetry(u);
                unit_features.T2Pdelay = t2pdelay(u);
                unit_features.T2Pratio = t2pratio(u);
                feature_names = fieldnames(unit_features);
                for f = 1:length(feature_names) %If slewrates could not be calculated // should happen very rarely
                   if isempty(unit_features.(feature_names{f}))
                       unit_features.(feature_names{f}) = 0;
                   end
                end
                waveform_features{u} = unit_features;
            end
        end
