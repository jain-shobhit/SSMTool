function prob = alg_construct_eqn(fhan, x0, pnames, p0)
%ALG_CONSTRUCT_EQN   Construct an instance of 'alg'.
%
% Add zero functions encoded in terms of problem variables and problem
% parameters, monitor functions that evaluate to the problem parameters,
% and corresponding inactive continuation parameters.
%
% PROB = ALG_CONSTRUCT_EQN(FHAN, X0, PNAMES, P0)
%
% PROB   - Continuation problem structure.
% FHAN   - Function handle to zero function.
% X0     - Initial solution guess for array of problem variables.
% PNAMES - Single string label or cell array of string labels for
%          continuation parameters tracking the value of the problem
%          parameters.
% P0     - Initial solution guess for array of problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_construct_eqn.m 2839 2015-03-05 17:09:01Z fschild $

xdim       = numel(x0);      % Number of problem variables
pdim       = numel(p0);      % Number of problem parameters
data.fhan  = fhan;           % Function handle to zero function
data.x_idx = (1:xdim)';      % Index set for problem variables
data.p_idx = xdim+(1:pdim)'; % Index set for problem parameters

prob = coco_prob();
prob = coco_add_func(prob, 'alg', @alg_F, data, 'zero', ...
  'u0', [x0(:); p0(:)]);
prob = coco_add_pars(prob, 'pars', data.p_idx, pnames);

end

function [data y] = alg_F(prob, data, u)
%ALG_F   COCO-compatible zero function wrapper.
%
% Supports zero functions encoded to expect arrays of problem variables
% and problem parameters, respectively, as input arguments.

x = u(data.x_idx); % Extract problem variables
p = u(data.p_idx); % Extract problem parameters

y = data.fhan(x, p);

end
