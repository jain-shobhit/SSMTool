function prob = po_close_orb(prob, tbid, data)
%PO_CLOSE_ORB   Append an instance of 'po' to problem.
%
% Add boundary and integral phase conditions to an instance of 'coll'.
%
% Differs from po_v1 by supporting the further subdivision of toolbox data
% implemented in coll_v4 and later, as well as adaptive remeshing.
%
% PROB = PO_CLOSE_ORB(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_close_orb.m 2839 2015-03-05 17:09:01Z fschild $

data.tbid = tbid;
data = coco_func_data(data); % Convert to func_data class for shared access
prob = coco_add_slot(prob, tbid, @po_update, data, 'update'); % Phase condition update function
segtbid      = coco_get_id(tbid, 'seg.coll'); % Create 'coll' toolbox instance identifier
[fdata uidx] = coco_get_func_data(prob, segtbid, 'data', 'uidx'); % Extract 'coll' data structure and context-dependent index array
prob = coco_add_func(prob, tbid, @po_F, @po_DFDU, data, 'zero', ...
  'uidx', uidx(fdata.maps.xbp_idx), 'remesh', @po_remesh);
fid  = coco_get_id(tbid, 'period');
prob = coco_add_pars(prob, fid, uidx(fdata.maps.T_idx), fid, 'active'); % Monitor period
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data y] = po_F(prob, data, u)
%PO_F   Evaluate boundary conditions and phase condition.
% 
% Integral phase condition is expressed as discretized linear operator.

x0 = u(data.x0_idx); % Extract trajectory end point at t=0
x1 = u(data.x1_idx); % Extract trajectory end point at t=1

y = [x0-x1; data.xp0*u];

end

function [data J] = po_DFDU(prob, data, u)
%PO_DFDU   Evaluate linearization of boundary conditions and phase condition.
% 
% Integral phase condition is expressed as discretized linear operator.

J = data.J;
end

function data = po_update(prob, data, cseg, varargin)
%PO_UPDATE   Update discretized linear operator.
%
% Use information about current solution to update discretization of linear
% operator corresponding to integral phase condition.

fid           = coco_get_id(data.tbid, 'seg.coll');
[fdata uidx]  = coco_get_func_data(prob, fid, 'data', 'uidx'); % 'coll' toolbox data and context-dependent index set
u             = cseg.src_chart.x; % Current chart
data.xp0      = u(uidx(fdata.maps.xbp_idx))'*data.intfac;
data.J(end,:) = data.xp0;

end

function [prob status xtr] = po_remesh(prob, data, chart, ub, vb)
%COLL_ERR_REMESH   Update dependency index set after remeshing.
%
% Associate a remesh action with the po_F function.

xtr    = []; % No invariant indices
data   = po_init_data(prob, data.tbid, data);      % Rebuild toolbox data
fid          = coco_get_id(data.tbid, 'seg.coll'); % Create 'coll' toolbox instance identifier
[fdata uidx] = coco_get_func_data(prob, fid, 'data', 'uidx'); % Extract 'coll' data and index array
prob   = coco_change_func(prob, data, 'uidx', uidx(fdata.maps.xbp_idx)); % Update data and solution
status = 'success';

end
