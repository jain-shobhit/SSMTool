function data = po_init_data(prob, tbid, data)
%PO_INIT_DATA   Initialize toolbox data for an instance of 'po'.
%
% Populate remaining fields of the toolbox data structure used by 'po'
% function objects.
%
% Differs from po_v1 by supporting the further subdivision of toolbox data
% implemented in coll_v4 and later.
%
% DATA = PO_INIT_DATA(PROB, TBID, DATA)
%
% DATA - Toolbox data structure.
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_init_data.m 2839 2015-03-05 17:09:01Z fschild $

stbid      = coco_get_id(tbid, 'seg.coll'); % Construct 'coll' toolbox instance identifier
[fdata u0] = coco_get_func_data(prob, stbid, 'data', 'u0'); % Extract 'coll' toolbox data and initial solution guess

dim  = fdata.int.dim;   % State-space dimension
NTST = fdata.coll.NTST; % Number of mesh intervals
NCOL = fdata.coll.NCOL; % Degree of interpolating polynomials
rows = [1:dim 1:dim];
cols = [fdata.maps.x0_idx fdata.maps.x1_idx];
vals = [ones(1,dim) -ones(1,dim)];
J    = sparse(rows, cols, vals, dim, dim*NTST*(NCOL+1)); % Jacobian of periodicity condition

data.x0_idx = fdata.maps.x0_idx;                   % Index array for trajectory end point at t=0
data.x1_idx = fdata.maps.x1_idx;                   % Index array for trajectory end point at t=1
data.intfac = fdata.maps.Wp'*fdata.mesh.wts2*fdata.maps.W;
data.xp0    = u0(fdata.maps.xbp_idx)'*data.intfac; % Discretized integral operator for phase conditions
data.J      = [J; data.xp0];                       % Linearization of boundary and phase conditions

end
