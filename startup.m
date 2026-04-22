% DeePhys startup script — adds all required paths.
%
% Usage: place this file on the MATLAB path, or call it explicitly from a script:
%   run('/path/to/DeePhys/startup.m')
%
% Alternatively, set the DEEPHYS_ROOT environment variable and call:
%   addpath(getenv('DEEPHYS_ROOT')); startup;

deephys_root = fileparts(mfilename('fullpath'));
addpath(genpath(deephys_root));
fprintf('DeePhys loaded from: %s\n', deephys_root);
