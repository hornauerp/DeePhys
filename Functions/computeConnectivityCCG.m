function connectivity = computeConnectivityCCG(spike_times, spike_units, n_units, sampling_rate, ccg_params)
% COMPUTECONNECTIVITYCCG  Compute CCG-based monosynaptic connectivity.
%
%   connectivity = computeConnectivityCCG(spike_times, spike_units, n_units, sampling_rate, ccg_params)
%
% INPUTS:
%   spike_times   - (N x 1) double, sorted spike times in seconds
%   spike_units   - (N x 1) double, unit assignment per spike
%   n_units       - (1 x 1) double, number of units
%   sampling_rate - (1 x 1) double, sampling rate in Hz
%   ccg_params    - struct: .BinSize, .Duration, .Conv_w, .Alpha
%
% OUTPUTS:
%   connectivity - struct with fields:
%       .ExcitatoryConnection  (E x 2) matrix of [source, target] pairs
%       .InhibitoryConnection  (I x 2) matrix
%       .bd                    (n_units x n_units) binary directed connectivity matrix
%       .CCGs                  (n_bins x n_units x n_units) cross-correlogram array

    arguments
        spike_times   (:,1) double
        spike_units   (:,1) double
        n_units       (1,1) double
        sampling_rate (1,1) double
        ccg_params    struct
    end

    [ccgR, ~] = CCG(spike_times, spike_units, ...
        'binSize', ccg_params.BinSize, ...
        'duration', ccg_params.Duration, ...
        'Fs', 1/sampling_rate);

    nCel = size(ccgR, 2);
    if nCel < n_units
        padded_ccg = zeros(size(ccgR,1), n_units, n_units);
        padded_ccg(:, 1:nCel, 1:nCel) = ccgR;
        ccgR = padded_ccg;
        nCel = n_units;
    end

    binSize = ccg_params.BinSize;
    conv_w  = ccg_params.Conv_w;
    alpha   = ccg_params.Alpha;

    nBins   = size(ccgR, 1);
    halfBin = round(nBins / 2);

    bonf_window = 0.005;
    pre_start   = 0.0032;
    post_start  = 0.0008;
    post_end    = 0.004;

    nBonf    = round(bonf_window / binSize) * 2;
    prebins  = round(halfBin - pre_start / binSize):halfBin;
    postbins = round(halfBin + post_start / binSize):round(halfBin + post_end / binSize);

    Pval   = nan(nBins, nCel, nCel);
    Pred   = zeros(nBins, nCel, nCel);
    Bounds = zeros(nBins, nCel, nCel, 2);

    n_pairs = nCel * (nCel - 1) / 2;
    sig_con     = zeros(n_pairs * 2, 2);
    sig_con_inh = zeros(n_pairs * 2, 2);
    n_exc = 0;
    n_inh = 0;

    for refcellID = 1:nCel
        for cell2ID = (refcellID+1):nCel
            cch = ccgR(:, refcellID, cell2ID);
            cchud = flipud(cch);

            [pvals, pred, ~] = bz_cch_conv(cch, conv_w);

            Pred(:, refcellID, cell2ID) = pred;
            Pval(:, refcellID, cell2ID) = pvals(:);
            Pred(:, cell2ID, refcellID) = flipud(pred(:));
            Pval(:, cell2ID, refcellID) = flipud(pvals(:));

            hiBound = poissinv(1 - alpha / nBonf, pred);
            loBound = poissinv(alpha / nBonf, pred);
            Bounds(:, refcellID, cell2ID, 1) = hiBound;
            Bounds(:, refcellID, cell2ID, 2) = loBound;
            Bounds(:, cell2ID, refcellID, 1) = flipud(hiBound(:));
            Bounds(:, cell2ID, refcellID, 2) = flipud(loBound(:));

            sig      = cch > hiBound;
            sig_inh  = cch < loBound;
            sigud     = flipud(sig);
            sigud_inh = flipud(sig_inh);

            sigpost = max(cch(postbins)) > poissinv(1 - alpha, max(cch(prebins)));
            if any(sig(postbins)) && sigpost
                n_exc = n_exc + 1;
                sig_con(n_exc, :) = [refcellID, cell2ID];
            end
            sigpre = max(cchud(postbins)) > poissinv(1 - alpha, max(cchud(prebins)));
            if any(sigud(postbins)) && sigpre
                n_exc = n_exc + 1;
                sig_con(n_exc, :) = [cell2ID, refcellID];
            end

            sigpost_inh = min(cch(postbins)) < poissinv(alpha, max(cch(prebins)));
            if any(sig_inh(postbins)) && sigpost_inh
                n_inh = n_inh + 1;
                sig_con_inh(n_inh, :) = [refcellID, cell2ID];
            end
            sigpre_inh = min(cchud(postbins)) < poissinv(alpha, max(cchud(prebins)));
            if any(sigud_inh(postbins)) && sigpre_inh
                n_inh = n_inh + 1;
                sig_con_inh(n_inh, :) = [cell2ID, refcellID];
            end
        end
    end

    sig_con     = sig_con(1:n_exc, :);
    sig_con_inh = sig_con_inh(1:n_inh, :);

    con_mat = zeros(nCel, nCel);
    for e = 1:size(sig_con, 1)
        con_mat(sig_con(e,1), sig_con(e,2)) = 1;
    end
    for i = 1:size(sig_con_inh, 1)
        con_mat(sig_con_inh(i,1), sig_con_inh(i,2)) = -1;
    end

    fprintf('Found %i excitatory and %i inhibitory connections\n', size(sig_con,1), size(sig_con_inh,1))

    connectivity.ExcitatoryConnection = sig_con;
    connectivity.InhibitoryConnection = sig_con_inh;
    connectivity.bd = con_mat;
    connectivity.CCGs = ccgR;
end
