function prob = alg_FO2FO(prob, oid, varargin)
%ALG_FO2FO   Append an instance of 'alg' and Moore-Spence system to problem.
%
% Support continuation of fold points from a previously located fold
% solution, stored to disk.
%
% PROB     = ALG_FO2FO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB     - Continuation problem structure.
% OID      - Target object instance identifier (string).
% RUN      - Run identifier (string).
% SOID     - Source object instance identifier (string).
% LAB      - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_FO2FO.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'alg');  % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run  = str.get;
if ischar(str.peek)
  soid = str.get;
else
  soid = oid;
end
lab = str.get;

[sol data] = alg_read_solution(soid, run, lab);  % Extract solution and toolbox data from disk
% The following function call is slightly different from Recipes for
% Continuation, 1st edition, page 433 ('off' instead of 'false').
prob       = coco_set(prob, tbid, 'FO', 'off');  % Turn off fold detection  
data       = alg_get_settings(prob, tbid, data); % Get toolbox settings
prob       = alg_construct_eqn(prob, tbid, data, sol); % Append 'alg' continuation problem

[data uidx] = coco_get_func_data(prob, tbid, 'data', 'uidx'); % Extract toolbox data structure and context-dependence index set
% The following function call is slightly different from Recipes for
% Continuation, 1st edition, page 433 (missing argument 'tbid').
prob        = alg_create_FO(prob, tbid, data, uidx, sol); % Append Moore-Spence system

end

% This function definition is slightly different from Recipes for
% Continuation, 1st edition, page 433 (missing argument 'tbid').
function prob = alg_create_FO(prob, tbid, data, uidx, sol)
%ALG_CREATE_FO   Append Moore-Spence system and initial nullvector to problem.
%
% PROB = ALG_CREATE_FO(PROB, TBID, DATA, UIDX, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - 'alg' toolbox data structure.
% UIDX - 'alg' context-dependent index array.
% SOL  - 'alg' solution.

Jx         = alg_fhan_DFDX(data, sol.x, sol.p);
[v0, ~]    = eigs(Jx, 1, 0); % Eigenvector corresponding to eigenvalue closest to 0
data.v_idx = data.p_idx(end)+(1:numel(v0)); % Index set for eigenvector

fid  = coco_get_id(tbid, 'fold_cond'); % Differs from Recipes for Continuation, 1st edition, page 433
prob = coco_add_func(prob, fid, @alg_FO, @alg_FO_DFDU, ...
  data, 'zero', 'uidx', [uidx(data.x_idx); uidx(data.p_idx)], ...
  'u0', v0);

end

function [data y] = alg_FO(prob, data, u)
%ALG_FO   Moore-Spence zero function.
%
% v is unit nullvector of Jacobian with respect to problem variables.

Jx = alg_fhan_DFDX(data, u(data.x_idx), u(data.p_idx));
v  = u(data.v_idx);
y  = [Jx*v; v'*v-1];

end

function [data J] = alg_FO_DFDU(prob, data, u)
%ALG_FO_DFDU   Linearization of Moore-Spence zero function.
%
% Second-order Jacobians of the zero function are here approximated by
% a mid-point finite-difference scheme applied to the first-order
% directional derivatives on the right-hand sides below:
%
% phi_{xx}(x,p)[t,.]=d/dh phi_{x}(x+h*t,p) at h=0
% phi_{xp}(x,p)[t,.]=d/dh phi_{p}(x+h*t,p) at h=0

x  = u(data.x_idx);
p  = u(data.p_idx);
Jx = alg_fhan_DFDX(data, x, p);
v  = u(data.v_idx);

h  = 1.0e-4*(1+norm(x));
J0 = alg_fhan_DFDX(data, x-h*v, p);
J1 = alg_fhan_DFDX(data, x+h*v, p);
Jxx = (0.5/h)*(J1-J0);

J0 = alg_fhan_DFDP(data, x-h*v, p);
J1 = alg_fhan_DFDP(data, x+h*v, p);
Jpx = (0.5/h)*(J1-J0);

J = [Jxx Jpx Jx; zeros(1,numel(data.x_idx)+numel(data.p_idx)) 2*v'];

end
