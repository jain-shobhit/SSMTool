function compalg_arg_check(prob, tbid, data)
%COMPALG_ARG_CHECK   Basic argument checking for 'compalg' toolbox.
%
% Validate user-supplied inputs and terminate execution with suitable error
% messages if the inputs fail to be of the correct type.
%
% Differs from compalg_v1 by validating against validated 'alg' instance
% data.
%
% COMPALG_ARG_CHECK(PROB, TBID, DATA)
%
% PROB     - continuation problem structure.
% TBID     - toolbox instance identifier.
% DATA     - initial toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: compalg_arg_check.m 2839 2015-03-05 17:09:01Z fschild $

assert(data.neqs~=0, '%s: insufficient number of equations', tbid);
pnum = [];
for i=1:data.neqs
  fid   = coco_get_id(tbid,sprintf('eqn%d.alg', i));
  fdata = coco_get_func_data(prob, fid, 'data');
  assert(isempty(fdata.pnames), ...
    '%s: parameter labels must not be passed to alg', tbid);
  assert(isempty(pnum) || pnum==numel(fdata.p_idx), '%s: %s', ...
    tbid, 'number of parameters must be equal for all equations');
  pnum = numel(fdata.p_idx);
end
assert(pnum==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''pnames''', ...
  tbid);

end
