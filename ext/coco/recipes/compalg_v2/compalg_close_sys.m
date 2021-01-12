function prob = compalg_close_sys(prob, tbid, data)
%COMPALG_CLOSE_SYS   Add gluing conditions to embedded instances of 'alg'.
%
% Extract context-dependent index arrays for the problem parameters
% associated with each instance of 'alg' and glue these together.
%
% Identical to compalg_v1.
%
% PROB = COMPALG_CLOSE_SYS(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: compalg_close_sys.m 2839 2015-03-05 17:09:01Z fschild $

sidx = cell(1, data.neqs);
for i=1:data.neqs % Context-dependent index arrays for copies of problem parameters
  stbid   = coco_get_id(tbid, sprintf('eqn%d.alg', i));
  [fdata uidx] = coco_get_func_data(prob, stbid, 'data', 'uidx');
  sidx{i} = uidx(fdata.p_idx);
end
for i=2:data.neqs % Glue redundant copies of problem parameters
  fid  = coco_get_id(tbid, sprintf('shared%d', i-1));
  prob = coco_add_glue(prob, fid, sidx{1}, sidx{i});
end
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  prob = coco_add_pars(prob, fid, sidx{1}, data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end
