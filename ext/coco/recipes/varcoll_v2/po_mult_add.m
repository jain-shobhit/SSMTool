function prob = po_mult_add(prob, segoid)
%PO_MULT_ADD   Add slot function computing Floquet multipliers of single-segment periodic orbit.
%
% Append slot function to a 'po' instance.
%
% Identical to varcoll_v1.
%
% PROB = PO_MULT_ADD(PROB, SEGOID)
%
% PROB   - Continuation problem structure.
% SEGOID - Segment object instance identifier.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_mult_add.m 2839 2015-03-05 17:09:01Z fschild $

data.tbid   = coco_get_id(segoid, 'var');
data.mnames = coco_get_id(segoid, 'multipliers'); % Column header
prob = coco_add_slot(prob, data.tbid, @po_mult_eigs_bddat, data, ...
  'bddat');

end

function [data res] = po_mult_eigs_bddat(prob, data, command, varargin)
%PO_MULT_EIGS_BDDAT   Slot function: add Floquet multipliers to bifurcation data.
%
% Extract Jacobian of time-T flow and compute its eigenvalues.
%
% Differs from varcoll_v1 by extracting the fundamental solution from the
% array of continuation variables.

switch command
  case 'init'
    res = {data.mnames};
  case 'data'
    [fdata uidx] = coco_get_func_data(prob, data.tbid, 'data', 'uidx'); % 'varcoll' data and index array
    chart = varargin{1}; % Current chart
    u     = chart.x(uidx);
    M     = reshape(u(fdata.ubp_idx), fdata.u_shp); % Fundamental solution
    M1    = M(fdata.M1_idx,:); % Jacobian of time-T flow
    res   = {eig(M1)}; % Includes trivial eigenvalue at 1
end

end
