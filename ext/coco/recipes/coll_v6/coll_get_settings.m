function data = coll_get_settings(prob, tbid, data)
%COLL_GET_SETTINGS   Read 'coll' toolbox instance settings.
%
% Merge user-supplied toolbox settings with default values.
%
% Differs from coll_v5 by support for varying discretization order.
%
% DATA = COLL_GET_SETTINGS(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data strcture.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_get_settings.m 2839 2015-03-05 17:09:01Z fschild $

defaults.NTST = 10; % Initial number of mesh intervals
defaults.NCOL = 4;  % Degree of interpolating polynomials
defaults.SAD  = 0.95; % Balance between error equidistribution and uniform mesh
if ~isfield(data, 'coll')
  data.coll = [];
end
data.coll = coco_merge(defaults, coco_merge(data.coll, ...
  coco_get(prob, tbid))); % Defaults < Stored < User-supplied
if ~coco_exist('TOL', 'class_prop', prob, tbid, '-no-inherit-all')
  data.coll.TOL = coco_get(prob, 'corr', 'TOL')^(2/3);
end % Guard against propagation of truncation errors
NTST = data.coll.NTST;
assert(numel(NTST)==1 && isnumeric(NTST) && mod(NTST,1)==0, ...
  '%s: input for option ''NTST'' is not an integer', tbid);
NCOL = data.coll.NCOL;
assert(numel(NCOL)==1 && isnumeric(NCOL) && mod(NCOL,1)==0, ...
  '%s: input for option ''NCOL'' is not an integer', tbid);

defaults.TOLINC = data.coll.TOL/5;          % Upper bound on adaptation window
defaults.TOLDEC = data.coll.TOL/20;         % Lower bound on adaptation window 
defaults.NTSTMN = min(  5, data.coll.NTST); % Lower bound on number of discretization intervals
defaults.NTSTMX = max(100, data.coll.NTST); % Upper bound on number of discretization intervals
data.coll = coco_merge(defaults, data.coll);

end
