function data = hspo_get_settings(prob, tbid)
%HSPO_GET_SETTINGS   Read 'hspo' toolbox instance settings.
%
% Merge user-supplied toolbox settings with default values.
%
% DATA = HSPO_GET_SETTINGS(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: hspo_get_settings.m 2839 2015-03-05 17:09:01Z fschild $

defaults.bifus = true;  % Detect saddle-node, period-doubling, and Neimark-Sacker bifurcations
defaults.NSad  = false; % Do not detect neutral saddles
data.hspo = coco_merge(defaults, coco_get(prob, tbid));
assert(islogical(data.hspo.bifus), ...
  '%s: input for ''bifus'' option is not boolean', tbid);
assert(islogical(data.hspo.NSad), ...
  '%s: input for ''NSad'' option is not boolean', tbid);

end
