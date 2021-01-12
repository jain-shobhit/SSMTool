function data = dft_get_settings(prob, tbid, data)
%DFT_GET_SETTINGS   Read 'dft' toolbox instance settings.
%
% Merge user-supplied toolbox settings with default values.
%
% DATA = DFT_GET_SETTINGS(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data strcture.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: dft_get_settings.m 2839 2015-03-05 17:09:01Z fschild $

defaults.NMAX = 8; % Maximum number of Fourier modes
defaults.NMIN = 3; % Minimum number of Fourier modes
defaults.NMOD = 3; % Initial number of Fourier modes
if ~isfield(data, 'dft')
  data.dft = [];
end
data.dft = coco_merge(defaults, coco_merge(data.dft, ...
  coco_get(prob, tbid))); % Defaults < Stored < User-supplied

if ~coco_exist('TOL', 'class_prop', prob, tbid, '-no-inherit-all')
  data.dft.TOL = coco_get(prob, 'corr', 'TOL')^(2/3);
end % Guard against propagation of truncation errors
defaults.TOLINC = data.dft.TOL/5;  % Upper bound on adaptation window
defaults.TOLDEC = data.dft.TOL/20; % Lower bound on adaptation window
data.dft = coco_merge(defaults, data.dft);

NMOD = data.dft.NMOD;
assert(numel(NMOD)==1 && isnumeric(NMOD) && mod(NMOD,1)==0, ...
  '%s: input for option ''NMOD'' is not an integer', tbid);
NMAX = data.dft.NMAX;
assert(numel(NMAX)==1 && isnumeric(NMAX) && mod(NMAX,1)==0, ...
  '%s: input for option ''NMAX'' is not an integer', tbid);
NMIN = data.dft.NMIN;
assert(numel(NMIN)==1 && isnumeric(NMIN) && mod(NMIN,1)==0, ...
  '%s: input for option ''NMIN'' is not an integer', tbid);
assert(NMIN<=NMOD && NMOD<=NMAX, ...
  '%s: input violates ''NMIN<=NMOD<=NMAX''', tbid);
end
