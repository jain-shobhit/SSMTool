function data = po_get_settings(prob, tbid, data)
%PO_GET_SETTINGS   Read 'po' toolbox instance settings.
%
% Merge user-supplied toolbox settings with default values.
%
% DATA = PO_GET_SETTINGS(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_get_settings.m 2839 2015-03-05 17:09:01Z fschild $

defaults.bifus = true;  % Detect saddle-node, period-doubling, and Neimark-Sacker bifurcations
defaults.NSad  = false; % Do not detect neutral saddles
data.po = coco_merge(defaults, coco_get(prob, tbid));
assert(islogical(data.po.bifus), ...
  '%s: input for ''bifus'' option is not boolean', tbid);
assert(islogical(data.po.NSad), ...
  '%s: input for ''NSad'' option is not boolean', tbid);

end
