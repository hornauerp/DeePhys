        function waveform_features = inferWaveformFeatures(obj, max_amplitudes, norm_wf_matrix)
            interpolation_factor = 10;
            ms_conversion = 10 * interpolation_factor; % 10 samples/ms * interpolation_factor = samples per ms in the interpolated waveform
            x = 1:size(norm_wf_matrix,1);
            xq = 1:(1/interpolation_factor):size(norm_wf_matrix,1);
            interp_wf_matrix = double(interp1(x,norm_wf_matrix,xq,'pchip'));
            zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %Function to detect zero crossings
            n_units = size(interp_wf_matrix,2);

            % Compute per-unit zero crossings via column-wise detection
            unit_zero_crossings = cell(1, n_units);
            for u = 1:n_units
                zc_idx = zci(interp_wf_matrix(:,u));
                unit_zero_crossings{u} = zc_idx(:)';
            end

            [unit_trough_value,unit_trough_idx] = min(interp_wf_matrix);
            waveform_features = cell(1,length(max_amplitudes));

            for u = 1:n_units
                unit_features = struct();

                % Per-unit peak cutouts using this unit's trough location
                p1_start = max(1, unit_trough_idx(u) - ms_conversion);
                p1_end   = unit_trough_idx(u);
                p2_start = unit_trough_idx(u);
                p2_end   = min(size(interp_wf_matrix,1), unit_trough_idx(u) + ms_conversion);

                peak_1_cutout = interp_wf_matrix(p1_start:p1_end, u);
                peak_2_cutout = interp_wf_matrix(p2_start:p2_end, u);

                [peak_1_value, ~] = max(peak_1_cutout);
                [peak_2_value, peak_2_rel_idx] = max(peak_2_cutout);
                peak_2_abs_idx = peak_2_rel_idx + unit_trough_idx(u) - 1;

                asymmetry_u = (peak_2_value - peak_1_value) / (peak_2_value + peak_1_value);
                t2pdelay_u  = (peak_2_abs_idx - unit_trough_idx(u)) / (obj.RecordingInfo.SamplingRate / 1000); %in [ms]
                t2pratio_u  = abs(unit_trough_value(u) / peak_2_value);

                % AUC features using per-unit zero crossings and waveform column
                zc = [1 unit_zero_crossings{u} length(xq)];
                [~,zc_pre_trough] = max(zc(zc<unit_trough_idx(u))); %Find zero crossing towards the trough
                zc_pre_trough = max([zc_pre_trough 2]);
                unit_features.AUC_peak_1 = trapz(interp_wf_matrix(zc(zc_pre_trough - 1):zc(zc_pre_trough), u));
                unit_features.AUC_trough = abs(trapz(interp_wf_matrix(zc(zc_pre_trough):zc(zc_pre_trough + 1), u)));
                if length(zc(zc_pre_trough:end)) > 3 %Ensure that there is another zero crossing at the end
                    unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):zc(zc_pre_trough + 2), u));
                else %Otherwise we integrate until the end of the waveform
                    unit_features.AUC_peak_2 = trapz(interp_wf_matrix(zc(zc_pre_trough + 1):size(interp_wf_matrix,1), u));
                end
                %We pad the signal with min/max values to ensure reliable slewrate
                %calculations
                rise_cutout = interp_wf_matrix(unit_trough_idx(u):peak_2_abs_idx,u);
                padded_rise = [ones(100,1)*rise_cutout(1); rise_cutout; ones(100,1)*rise_cutout(end)];
                unit_features.Rise = mean(slewrate(padded_rise,interpolation_factor));
                decay_cutout = interp_wf_matrix(peak_2_abs_idx:end,u);
                decay_cutout = decay_cutout(1:find(decay_cutout==min(decay_cutout), 1, 'first'));
                padded_decay = [ones(100,1)*decay_cutout(1); decay_cutout; ones(100,1)*decay_cutout(end)];
                unit_features.Decay = mean(slewrate(padded_decay,interpolation_factor));
                unit_features.Asymmetry = asymmetry_u;
                unit_features.T2Pdelay = t2pdelay_u;
                unit_features.T2Pratio = t2pratio_u;
                feature_names = fieldnames(unit_features);
                for f = 1:length(feature_names) %If slewrates could not be calculated // should happen very rarely
                   if isempty(unit_features.(feature_names{f}))
                       unit_features.(feature_names{f}) = 0;
                   end
                end
                waveform_features{u} = unit_features;
            end
        end
