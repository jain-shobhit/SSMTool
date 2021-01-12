function prob = po_mult_add(prob, segoid)
%PO_MULT_ADD   Add monitor function computing Floquet multipliers of single-segment periodic orbit.
%
% Append monitor function to a 'po' instance.
%
% PROB = PO_MULT_ADD(PROB, SEGOID)
%
% PROB   - Continuation problem structure.
% SEGOID - Segment object instance identifier.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_mult_add.m 2871 2015-08-03 12:34:16Z hdankowicz $

data.coid   = coco_get_id(segoid, 'coll');
data.tbid   = coco_get_id(segoid, 'coll.test');
data.mnames = coco_get_id(segoid, 'multipliers'); % Column header 
fdata = coco_get_func_data(prob, data.coid, 'data');
mfid  = coco_get_def_par_names(data.mnames, 1:fdata.coll_seg.int.dim);
prob  = coco_add_func(prob, data.mnames, @floquet, data, ...
  'regular', mfid, 'requires', data.tbid, 'remesh', @remesh);

end

function [data, y] = floquet(prob, data, u) %#ok<INUSD>
%FLOQUET   Compute Floquet multipliers.
%
% Extract Jacobian of time-T flow and compute its eigenvalues.

fdata = coco_get_func_data(prob, data.tbid, 'data');
ctst  = fdata.coll_tst;
M0    = ctst.M(ctst.M0_idx,:);
M1    = ctst.M(ctst.M1_idx,:);
y     = eig(full(M1/M0)); % Includes trivial eigenvalue at 1
    
end

function [prob, stat, xtr] = remesh(prob, data, chart, ub, Vb) %#ok<INUSD>

xtr  = [];
prob = coco_change_func(prob, data);
stat = 'success';

end
