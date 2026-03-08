       function inferConnectivityCCG(obj,ccg_type)
           arguments
            obj MEArecording
            ccg_type string = "CCG" %or "FullCCG"
           end
           fh = str2func("calculate" + ccg_type);
           [ccgR1,tR] = fh(obj);
           if size(ccgR1,2) < length(obj.Units) %Pad if the last unit(s) have no spikes
               padded_ccg = zeros(size(ccgR1,1),length(obj.Units),length(obj.Units));
               padded_ccg(:,1:size(ccgR1,2),1:size(ccgR1,2)) = ccgR1;
               ccgR1 = padded_ccg;
           end
           binSize = obj.Parameters.CCG.BinSize; %.5ms
           conv_w = obj.Parameters.CCG.Conv_w;  % 10ms window
           alpha = obj.Parameters.CCG.Alpha; %high frequency cut off, must be .001 for causal p-value matrix
%            Fs = 1/obj.RecordingInfo.SamplingRate;

           nCel=size(ccgR1,2);

           ccgR = nan(size(ccgR1,1),nCel,nCel);
           ccgR(:,1:size(ccgR1,2),1:size(ccgR1,2)) = ccgR1;

           % get  CI for each CCG
           Pval=nan(length(tR),nCel,nCel);
           Pred=zeros(length(tR),nCel,nCel);
           Bounds=zeros(size(ccgR,1),nCel,nCel);
           sig_con = [];
           sig_con_inh = [];

           for refcellID=1:max(nCel)
               for cell2ID=1:max(nCel)

                   if(refcellID==cell2ID)
                       continue;
                   end

                   cch=ccgR(:,refcellID,cell2ID);			% extract corresponding cross-correlation histogram vector

                   [pvals,pred,~]=bz_cch_conv(cch,conv_w);
                   % Store predicted values and pvalues for subsequent plotting
                   Pred(:,refcellID,cell2ID)=pred;
                   Pval(:,refcellID,cell2ID)=pvals(:);
                   Pred(:,cell2ID,refcellID)=flipud(pred(:));
                   Pval(:,cell2ID,refcellID)=flipud(pvals(:));

                   % Calculate upper and lower limits with bonferonni correction
                   % monosynaptic connection will be +/- 4 ms

                   nBonf = round(.005/binSize)*2;
                   hiBound=poissinv(1-alpha/nBonf,pred);
                   loBound=poissinv(alpha/nBonf, pred);
                   Bounds(:,refcellID,cell2ID,1)=hiBound;
                   Bounds(:,refcellID,cell2ID,2)=loBound;

                   Bounds(:,cell2ID,refcellID,1)=flipud(hiBound(:));
                   Bounds(:,cell2ID,refcellID,2)=flipud(loBound(:));

                   sig = cch>hiBound;
                   sig_inh= cch < loBound;

                   % Find if significant periods falls in monosynaptic window +/- 4ms
                   prebins = round(length(cch)/2 - .0032/binSize):round(length(cch)/2);
                   postbins = round(length(cch)/2 + .0008/binSize):round(length(cch)/2 + .004/binSize);
                   cchud  = flipud(cch);
                   sigud  = flipud(sig);
                   sigud_inh=flipud(sig_inh);

                   sigpost=max(cch(postbins))>poissinv(1-alpha,max(cch(prebins)));
                   sigpre=max(cchud(postbins))>poissinv(1-alpha,max(cchud(prebins)));

                   sigpost_inh=min(cch(postbins))<poissinv(alpha,mean(cch(prebins)));
                   sigpre_inh=min(cchud(postbins))<poissinv(alpha,mean(cchud(prebins)));
                   %check which is bigger
                   if (any(sigud(prebins)) && sigpre)
                       %test if causal is bigger than anti causal
                       sig_con = [sig_con;cell2ID refcellID];
                   end

                   if (any(sig(postbins)) && sigpost)
                       sig_con = [sig_con;refcellID cell2ID];
                   end

                   if (any(sigud_inh(prebins)) && sigpre_inh)
                       %test if causal is bigger than anti causal
                       sig_con_inh = [sig_con_inh;cell2ID refcellID];
                   end
                   if (any(sig_inh(postbins)) && sigpost_inh)
                       sig_con_inh = [sig_con_inh;refcellID cell2ID];
                   end
               end

           end

           sig_con=unique(sig_con,'rows');
           sig_con_inh=unique(sig_con_inh, 'rows');

           if(~isempty(sig_con))
               ccg_vec=[];
               for jj=1:size(sig_con,1)
                   ccg_vec=[ccg_vec ccgR(:,sig_con(jj,1),sig_con(jj,2))];
               end
           else
               ccg_vec=[];
           end

           if(~isempty(sig_con_inh))
               ccg_vec_inh=[];
               for jj=1:size(sig_con_inh,1)
                   ccg_vec_inh=[ccg_vec_inh ccgR(:,sig_con_inh(jj,1),sig_con_inh(jj,2))];
               end
           else
               ccg_vec_inh=[];
           end
           fprintf('Found %i excitatory and %i inhibitory connections\n', size(sig_con,1), size(sig_con_inh,1))
           con_mat = zeros(size(ccgR,[2,3]));
           for e = 1:size(sig_con,1)
               con_mat(sig_con(e,1),sig_con(e,2)) = 1;
           end
           for i = 1:size(sig_con_inh,1)
               con_mat(sig_con_inh(i,1),sig_con_inh(i,2)) = -1;
           end

%            silent_units = find(arrayfun(@(x) isempty(x.SpikeTimes),obj.Units));

           obj.Connectivity.(ccg_type).ExcitatoryConnection = sig_con;
           obj.Connectivity.(ccg_type).InhibitoryConnection = sig_con_inh;
           obj.Connectivity.(ccg_type).bd = con_mat;
           obj.Connectivity.(ccg_type).CCGs = ccgR;
       end
