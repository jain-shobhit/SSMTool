function prob = po_mult_add(prob, segoid)
%PO_MULT_ADD   Add slot and/or monitor function computing Floquet multipliers of single-segment periodic orbit.
%
% Append slot and/or monitor function to a 'po' instance.
%
% PROB = PO_MULT_ADD(PROB, SEGOID)
%
% PROB   - Continuation problem structure.
% SEGOID - Segment object instance identifier.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_mult_add.m 2839 2015-03-05 17:09:01Z fschild $

data.coid   = coco_get_id(segoid, 'coll'); % Not in Recipes for Continuation, 1st edition.
data.tbid   = coco_get_id(segoid, 'var');

% This is the solution from Recipes for Continuation, 1st edition, page
% 258, albeit with the corrected implementation of po_mult_eigs_bddat,
% which eliminates the incorrect use of function data.

data.mnames = coco_get_id(segoid, 'multipliers'); % Column header 
prob = coco_add_slot(prob, data.tbid, @po_mult_eigs_bddat, data, ...
  'bddat');

% This is an alternative solution that makes correct use of function data.

fdata = coco_get_func_data(prob, data.coid, 'data');
mfid  = coco_get_def_par_names(data.mnames, 1:fdata.dim);
prob  = coco_add_func(prob, data.mnames, @po_mult_eigs, data, 'regular', ...
  mfid, 'requires', data.tbid); % Explicit dependence on 'varcoll'

end

function [data y] = po_mult_eigs(prob, data, u)
%PO_MULT_EIGS   Compute Floquet multipliers.
%
% This function is not included in Recipes for Continuation, 1st edition.
%
% Extract Jacobian of time-T flow and compute its eigenvalues.

fdata = coco_get_func_data(prob, data.tbid, 'data');
M1    = fdata.M(fdata.M1_idx,:);
y     = eig(full(M1)); % Includes trivial eigenvalue at 1
    
end

function [data res] = po_mult_eigs_bddat(prob, data, command, varargin)
%PO_MULT_EIGS_BDDAT   Slot function: add Floquet multipliers to bifurcation data.
%
% Extract Jacobian of time-T flow and compute its eigenvalues.

switch command
  case 'init'
    res   = {data.mnames};
  case 'data'
    % The following two lines from Recipes for Continuation, 1st edition,
    % page 257, make incorrect use of function data to extract the Jacobian
    % of the time-T flow. Notably, the available function data of the
    % instance of the 'varcoll' toolbox need not be associated with the
    % chart that is being flushed to the bifurcation data.
    
    % fdata = coco_get_func_data(prob, data.tbid, 'data');
    % M1    = fdata.M(fdata.M1_idx,:);
 
    % The following four lines provide a correct solution to the problem of
    % extracting the Jacobian of the time-T flow. The implementation
    % guarantees a just-in-time solution of the variational equation.
    
    [fdata uidx] = coco_get_func_data(prob, data.coid, 'data', 'uidx');
    chart = varargin{1};           	% Extract the current chart
    u     = chart.x(uidx);
    M1    = var_eval_sol(fdata, u); % Solve the variational equation
    
    res   = {eig(full(M1))};        % Includes trivial eigenvalue at 1
end

end
