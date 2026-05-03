        function Burst = detectBurstsISIN(obj, N, ISI_N)
            %Adapted from Bakkum et al., 2014
            if N == 0
                Burst.T_start = [];
                Burst.T_end = [];
            else
                tic;
                [spike_times, ~] = obj.removeTonic();
                dT = zeros(N,length(spike_times))+inf;
                for j = 0:N-1
                    dT(j+1,N:length(spike_times)-(N-1)) = spike_times( (N:end-(N-1))+j ) - ...
                        spike_times( (1:end-(N-1)*2)+j );
                end
                Criteria = zeros(size(spike_times)); % Initialize to zero
                Criteria( min(dT)<=ISI_N ) = 1; % Spike passes condition if it is
                % included in a set of N spikes
                % with ISI_N <= threshold.
                % %% Assign burst numbers to each spike
                SpikeBurstNumber = zeros(size(spike_times)) - 1; % Initialize to '-1'
                INBURST = 0; % In a burst (1) or not (0)
                NUM_ = 0; % Burst Number iterator
                NUMBER = -1; % Burst Number assigned
                BL = 0; % Burst Length
                for i = N:length(spike_times)
                    if INBURST == 0 % Was not in burst.
                        if Criteria(i) % Criteria met, now in new burst.
                            INBURST = 1; % Update.
                            NUM_ = NUM_ + 1;
                            NUMBER = NUM_;
                            BL = 1;
                        else % Still not in burst, continue.
                            continue
                        end
                    else % Was in burst.
                        if ~ Criteria(i) % Criteria no longer met.
                            INBURST = 0; % Update.
                            if BL<N % Erase if not big enough.
                                SpikeBurstNumber(SpikeBurstNumber==NUMBER) = -1;
                                NUM_ = NUM_ - 1;
                            end
                            NUMBER = -1;
                        elseif diff(spike_times([i-(N-1) i])) > ISI_N && BL >= N
                            % This conditional statement is necessary to split apart
                            % consecutive bursts that are not interspersed by a tonic spike
                            % (i.e. Criteria == 0). Occasionally in this case, the second
                            % burst has fewer than 'N' spikes and is therefore deleted in
                            % the above conditional statement (i.e. 'if BL<N').
                            %
                            % Skip this if at the start of a new burst (i.e. 'BL>=N'
                            % requirement).
                            %
                            NUM_ = NUM_ + 1; % New burst, update number.
                            NUMBER = NUM_;
                            BL = 1; % Reset burst length.
                        else % Criteria still met.
                            BL = BL + 1; % Update burst length.
                        end
                    end
                    SpikeBurstNumber(i) = NUMBER; % Assign a burst number to
                    % each spike.
                end
                % %% Assign Burst information
                MaxBurstNumber = max(SpikeBurstNumber);
                Burst.T_start = zeros(1,MaxBurstNumber); % Burst start time [sec]
                Burst.T_end = zeros(1,MaxBurstNumber); % Burst end time [sec]
                for i = 1:MaxBurstNumber
                    ID = find( SpikeBurstNumber==i );
                    Burst.T_start(i) = spike_times(ID(1));
                    Burst.T_end(i) = spike_times(ID(end));
                end
                run_time = toc;
                fprintf('Detected %i bursts in %.2f seconds using %0.2f minutes of spike data.\n', ...
                    MaxBurstNumber,run_time,diff(spike_times([1 end]))/60);
            end
        end
