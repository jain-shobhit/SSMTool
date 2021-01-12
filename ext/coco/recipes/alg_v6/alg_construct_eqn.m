function prob = alg_construct_eqn(prob, tbid, data, sol)
%ALG_CONSTRUCT_EQN   Append an instance of 'alg' to problem.
%
% Add zero functions encoded in terms of problem variables and problem
% parameters, monitor functions that evaluate to the problem parameters,
% and corresponding inactive continuation parameters.
%
% Differs from alg_v5 by the use of an optional slot function in response
% to the 'bddat' signal to add content to the bifurcation data cell array.
%
% PROB = ALG_CONSTRUCT_EQN(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_construct_eqn.m 2839 2015-03-05 17:09:01Z fschild $

prob = coco_add_func(prob, tbid, @alg_F, @alg_DFDU, data, 'zero', ...
  'u0', sol.u);
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
  prob = coco_add_pars(prob, fid, uidx(data.p_idx), data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');
if data.alg.norm % Optional slot function
  data.tbid = tbid;
  prob = coco_add_slot(prob, tbid, @alg_bddat, data, 'bddat');
end

end

function [data y] = alg_F(prob, data, u)
%ALG_F   COCO-compatible zero function wrapper.
%
% Supports zero functions encoded to expect two arrays of problem variables
% and problem parameters, respectively, as input arguments.
%
% Identical to alg_v1. 

x = u(data.x_idx); % Extract problem variables
p = u(data.p_idx); % Extract problem parameters

y = data.fhan(x, p);

end

function [data J] = alg_DFDU(prob, data, u)
%ALG_DFDU   COCO-compatible linearization of zero function wrapper.
%
% Supports linearization of zero functions encoded to expect two arrays of
% problem variables and problem parameters, respectively, as input
% arguments. When dfdxhan and/or dfdphan are empty, approximate Jacobians
% are obtained using numerical differentiation.
%
% Identical to alg_v2. 

x = u(data.x_idx); % Extract problem variables
p = u(data.p_idx); % Extract problem parameters

if isempty(data.dfdxhan)
  J1 = coco_ezDFDX('f(x,p)',data.fhan, x, p);
else
  J1 = data.dfdxhan(x, p);
end
if isempty(data.dfdphan)
  J2 = coco_ezDFDP('f(x,p)',data.fhan, x, p);
else
  J2 = data.dfdphan(x, p);
end
J = sparse([J1 J2]);

end

function [data res] = alg_bddat(prob, data, command, varargin)
%ALG_BDDAT   Append to bifurcation data cell array.
%
% Function adds column containing Euclidean norm of array of
% problem variables to bifurcation data cell array.

res = {};
switch command
  case 'init'
    % The following line corrects a missing toolbox instance identifier in
    % Recipes for Continuation, 1st edition, page 108.
    res   = sprintf('||%s.x||', data.tbid); % Column header
  case 'data'
    chart = varargin{1}; % Current chart
    uidx  = coco_get_func_data(prob, data.tbid, 'uidx');
    u     = chart.x(uidx);
    res   = norm(u(data.x_idx), 2); % Column data
end

end
