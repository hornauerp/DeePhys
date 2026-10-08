function tf = acgColumns(feat_names)
% ACGCOLUMNS  True for the ACG feature columns (FullACG* / ACG*); the rest are waveform.
    feat_names = string(feat_names);
    tf = startsWith(feat_names, "FullACG") | startsWith(feat_names, "ACG");
end
