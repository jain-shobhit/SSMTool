function data = alg_get_settings(prob, tbid, data)
%ALG_GET_SETTINGS   Read 'alg' toolbox instance settings.
%
% Merge user-supplied toolbox settings with default values.
%
% DATA = ALG_GET_SETTINGS(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
%
% See also: coco_set, coco_get, coco_merge

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_get_settings.m 2839 2015-03-05 17:09:01Z fschild $

defaults.norm = false; % Do not append Euclidean norm to bifurcation data
data.alg = coco_merge(defaults, coco_get(prob, tbid));
assert(islogical(data.alg.norm), ...
  '%s: input for ''norm'' option is not boolean', tbid);

end
