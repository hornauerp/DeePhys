       function inferConnectivityCCG(obj,ccg_type)
           arguments
            obj MEArecording
            ccg_type string = "CCG" %or "FullCCG"
           end
           % Compute CCG (single call — calculateCCG does not cache its results)
           fh = str2func("calculate" + ccg_type);
           [ccgR,tR] = fh(obj);

           nCel = size(ccgR, 2);
           if nCel < length(obj.Units)
               padded_ccg = zeros(size(ccgR,1), length(obj.Units), length(obj.Units));
               padded_ccg(:,1:nCel,1:nCel) = ccgR;
               ccgR = padded_ccg;
               nCel = length(obj.Units);
           end

           binSize = obj.Parameters.CCG.BinSize;
           conv_w = obj.Parameters.CCG.Conv_w;
           alpha = obj.Parameters.CCG.Alpha;

           nBins = size(ccgR, 1);
           halfBin = round(nBins / 2);

           % Monosynaptic window parameters (seconds)
           bonf_window = 0.005;        % Bonferroni correction window
           pre_start   = 0.0032;       % pre-synaptic window start offset
           post_start  = 0.0008;       % post-synaptic window start offset
           post_end    = 0.004;        % post-synaptic window end offset

           % Pre-compute bin indices for monosynaptic window
           nBonf    = round(bonf_window / binSize) * 2;
           prebins  = round(halfBin - pre_start / binSize):halfBin;
           postbins = round(halfBin + post_start / binSize):round(halfBin + post_end / binSize);

           % Pre-allocate output arrays
           Pval   = nan(nBins, nCel, nCel);
           Pred   = zeros(nBins, nCel, nCel);
           Bounds = zeros(nBins, nCel, nCel, 2);

           % Pre-allocate connection lists (upper bound: each pair produces at most 2 connections)
           n_pairs = nCel * (nCel - 1) / 2;
           sig_con     = zeros(n_pairs * 2, 2);
           sig_con_inh = zeros(n_pairs * 2, 2);
           n_exc = 0;
           n_inh = 0;

           % Iterate upper triangle only — each pair's CCG gives both causal directions
           for refcellID = 1:nCel
               for cell2ID = (refcellID+1):nCel
                   cch = ccgR(:, refcellID, cell2ID);
                   cchud = flipud(cch);

                   % Convolution-based prediction and p-values
                   [pvals, pred, ~] = bz_cch_conv(cch, conv_w);

                   % Store for both directions (symmetric via flip)
                   Pred(:, refcellID, cell2ID) = pred;
                   Pval(:, refcellID, cell2ID) = pvals(:);
                   Pred(:, cell2ID, refcellID) = flipud(pred(:));
                   Pval(:, cell2ID, refcellID) = flipud(pvals(:));

                   % Bonferroni-corrected bounds
                   hiBound = poissinv(1 - alpha / nBonf, pred);
                   loBound = poissinv(alpha / nBonf, pred);
                   Bounds(:, refcellID, cell2ID, 1) = hiBound;
                   Bounds(:, refcellID, cell2ID, 2) = loBound;
                   Bounds(:, cell2ID, refcellID, 1) = flipud(hiBound(:));
                   Bounds(:, cell2ID, refcellID, 2) = flipud(loBound(:));

                   % Significance masks
                   sig      = cch > hiBound;
                   sig_inh  = cch < loBound;
                   sigud     = flipud(sig);
                   sigud_inh = flipud(sig_inh);

                   % --- Excitatory connections ---
                   % Forward: ref → cell2 (post-synaptic peak in cch)
                   sigpost = max(cch(postbins)) > poissinv(1 - alpha, max(cch(prebins)));
                   if any(sig(postbins)) && sigpost
                       n_exc = n_exc + 1;
                       sig_con(n_exc, :) = [refcellID, cell2ID];
                   end
                   % Reverse: cell2 → ref (pre-synaptic peak = post-synaptic in flipped cch)
                   sigpre = max(cchud(postbins)) > poissinv(1 - alpha, max(cchud(prebins)));
                   if any(sigud(prebins)) && sigpre
                       n_exc = n_exc + 1;
                       sig_con(n_exc, :) = [cell2ID, refcellID];
                   end

                   % --- Inhibitory connections ---
                   sigpost_inh = min(cch(postbins)) < poissinv(alpha, mean(cch(prebins)));
                   if any(sig_inh(postbins)) && sigpost_inh
                       n_inh = n_inh + 1;
                       sig_con_inh(n_inh, :) = [refcellID, cell2ID];
                   end
                   sigpre_inh = min(cchud(postbins)) < poissinv(alpha, mean(cchud(prebins)));
                   if any(sigud_inh(prebins)) && sigpre_inh
                       n_inh = n_inh + 1;
                       sig_con_inh(n_inh, :) = [cell2ID, refcellID];
                   end
               end
           end

           % Trim pre-allocated arrays
           sig_con     = sig_con(1:n_exc, :);
           sig_con_inh = sig_con_inh(1:n_inh, :);

           % Extract CCG vectors for significant connections (vectorized)
           ccgR_2d = reshape(ccgR, nBins, []);  % nBins × (nCel*nCel)
           if ~isempty(sig_con)
               idx = sub2ind([nCel, nCel], sig_con(:,1), sig_con(:,2));
               ccg_vec = ccgR_2d(:, idx);
           else
               ccg_vec = [];
           end

           if ~isempty(sig_con_inh)
               idx_inh = sub2ind([nCel, nCel], sig_con_inh(:,1), sig_con_inh(:,2));
               ccg_vec_inh = ccgR_2d(:, idx_inh);
           else
               ccg_vec_inh = [];
           end

           fprintf('Found %i excitatory and %i inhibitory connections\n', size(sig_con,1), size(sig_con_inh,1))

           % Build binary connection matrix
           con_mat = zeros(nCel, nCel);
           for e = 1:size(sig_con, 1)
               con_mat(sig_con(e,1), sig_con(e,2)) = 1;
           end
           for i = 1:size(sig_con_inh, 1)
               con_mat(sig_con_inh(i,1), sig_con_inh(i,2)) = -1;
           end

           obj.Connectivity.(ccg_type).ExcitatoryConnection = sig_con;
           obj.Connectivity.(ccg_type).InhibitoryConnection = sig_con_inh;
           obj.Connectivity.(ccg_type).bd = con_mat;
           obj.Connectivity.(ccg_type).CCGs = ccgR;
       end
