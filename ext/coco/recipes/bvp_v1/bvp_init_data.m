function data = bvp_init_data(prob, tbid, data)
%BVP_INIT_DATA   Initialize toolbox data for an instance of 'bvp'.
%
% Populate remaining fields of the toolbox data structure used by 'bvp'
% function objects.
%
% DATA = BVP_INIT_DATA(PROB, TBID, DATA)
%
% DATA - Toolbox data structure.
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: bvp_init_data.m 2839 2015-03-05 17:09:01Z fschild $

segtbid = coco_get_id(tbid, 'seg.coll'); % Construct 'coll' toolbox instance identifier
fdata   = coco_get_func_data(prob, segtbid, 'data'); % Extract 'coll' toolbox data

data.T_idx  = 1;                              % Index for interval length
data.x0_idx = 1+(1:fdata.dim)';               % Index array for trajectory end point at t=0
data.x1_idx = 1+fdata.dim+(1:fdata.dim)';     % Index array for trajectory end point at t=1
data.p_idx  = 1+2*fdata.dim +(1:fdata.pdim)'; % Index array for problem parameters

end
