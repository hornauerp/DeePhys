function [oscillatory,fractal, original] = foof(data,ops)
arguments
    data = struct
    ops.freq_range (1,2) double = [2, 100]
end
cfg = [];
data = ft_preprocessing(cfg,data);
% b_lengths = cellfun(@(x) x(end)-x(1),data.trial);
% compute the fractal and original spectra
cfg               = [];
cfg.foilim        = ops.freq_range;
% cfg.foi = ops.freq_range(1):ops.freq_range(2);
cfg.pad           = ceil((data.sampleinfo(1,2) - data.sampleinfo(1,1))/data.fsample) + 1;
% if length(data.trial) == 1
    cfg.tapsmofrq = 1;%1000/min(cellfun(@length, data.trial));
% else
    % cfg.taper = 'hanning';
% end
cfg.method        = 'mtmfft';
cfg.output        = 'fooof_aperiodic';
fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'pow';
original = ft_freqanalysis(cfg, data);

cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2./x1';  % equivalent to 10^(log10(x2)-log10(x1))
oscillatory = ft_math(cfg, fractal, original);
end