       function [inburst_activity, inburst_units] =  getBurstStatistics(obj, ops)
           arguments
               obj MEArecording
               ops.binning = 0.01
           end
           if isempty(obj.Bursts) %Burst detection was not performed
               detectBursts(obj, merge_factor=obj.Parameters.Bursts.MergeFactor);
           end
           if length(obj.Bursts.T_start) < 3 %Not enough bursts were detected — use NaN, not 0
               Burst.MeanInterBurstInterval = NaN;
               Burst.BurstDuration = NaN;
               Burst.RiseTime = NaN;
               Burst.FallTime = NaN;
               Burst.PeakFR = NaN;
               Burst.StdFR = NaN;
               Burst.StdBurstDuration = NaN;
               Burst.StdInterBurstInterval = NaN;
               Burst.IntraFiringRate = NaN;
               Burst.InterFiringRate = NaN;
           else
               IBIs = obj.Bursts.T_start(2:end) - obj.Bursts.T_end(1:end-1);
               BDs = obj.Bursts.T_end - obj.Bursts.T_start;

               [RiseTime, FallTime, PeakFR, StdFR] = deal(nan(1,length(BDs)));
               Fs = round(1/ops.binning);
               inburst_activity = cell(1,length(BDs));
               inburst_units = cell(1,length(BDs));
               for iburst = 1:length(BDs)
                   try
                       burst_activity = obj.Spikes.Times(obj.Spikes.Times <= obj.Bursts.T_end(iburst) & obj.Spikes.Times >= obj.Bursts.T_start(iburst));
                       binned_activity = histcounts(burst_activity,'BinWidth',ops.binning);
                       [max_binned, imax_binned] = max(binned_activity);
                       norm_activity = smoothdata(binned_activity/max_binned,'lowess');
                       RiseTime(iburst) = imax_binned * ops.binning;
                       FallTime(iburst) = (length(binned_activity) - imax_binned) * ops.binning;
                       %%% Old approch, but seems too unreliable %%%

                       % if length(norm_activity) < 4
                       %     [peak_1,peak_1_idx] = max(norm_activity);
                       %     peak_2 = peak_1; peak_2_idx = peak_1_idx;
                       % else
                       %     [peak_1,peak_1_idx] = findpeaks(norm_activity,"MinPeakProminence",0.3,NPeaks=1);
                       %     [peak_2,peak_2_idx] = findpeaks(norm_activity(end:-1:1),"MinPeakProminence",0.3,NPeaks=1);
                       %     peak_2_idx = length(norm_activity) - peak_2_idx + 1;
                       % end
                       % if isempty(peak_1)
                       %     [peak_1,peak_1_idx] = max(norm_activity);
                       % end
                       % if isempty(peak_2)
                       %     [peak_2,peak_2_idx] = max(norm_activity);
                       % end
                       %
                       % pre_lower_level = min(norm_activity(1:peak_1_idx));
                       % post_lower_level = min(norm_activity(peak_2_idx:end));
                       % % Using mean to account for the rare case where rise/fall is split
                       % RiseTime(iburst) = mean(risetime([ones(1,Fs)*pre_lower_level norm_activity(1:peak_1_idx) ones(1,Fs)*peak_1],Fs));
                       % FallTime(iburst) = mean(falltime([ones(1,Fs)*peak_2 norm_activity(peak_2_idx:end) ones(1,Fs)*post_lower_level],Fs));

                       %%% Old approch, but seems too unreliable %%%

                       PeakFR(iburst) = max_binned./length(obj.Units) * Fs;
                       StdFR(iburst) = std(norm_activity);
                   catch ME
                       warning(ME.message)
                   end

               end
               [inburst, outburst] = obj.getBurstSpikes();
               inter_spikes = cellfun(@length, outburst.spikes);
               inter_fr = inter_spikes(2:end-1) ./ (IBIs * length(obj.Units));
               intra_spikes = cellfun(@length, inburst.spikes);
               intra_fr = intra_spikes ./ (BDs * length(obj.Units));

               cellfun_input = {IBIs, BDs, RiseTime, FallTime, PeakFR, StdFR};

               % Perform outlier removal if method is specified
               if ~isempty(obj.Parameters.Outlier.Method)
                   cellfun_input = cellfun(@(x) rmoutliers(x,obj.Parameters.Outlier.Method,'ThresholdFactor',obj.Parameters.Outlier.ThresholdFactor),...
                       cellfun_input,'un',0);
               end
               feature_means = cellfun(@(x) median(x,'omitnan'), cellfun_input,'un',0);

               [Burst.MeanInterBurstInterval,...
                   Burst.BurstDuration,...
                   Burst.RiseTime,...
                   Burst.FallTime,...
                   Burst.PeakFR,...
                   Burst.StdFR] = feature_means{:};
               Burst.StdBurstDuration = std(BDs);
               Burst.StdInterBurstInterval = std(IBIs);
               Burst.IntraFiringRate = median(intra_fr);
               Burst.InterFiringRate = median(inter_fr);
           end
           obj.NetworkFeatures.Burst = struct2table(Burst);
       end
