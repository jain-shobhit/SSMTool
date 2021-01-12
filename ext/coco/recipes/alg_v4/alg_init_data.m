function data = alg_init_data(data, x0, p0)
%ALG_INIT_DATA   Initialize toolbox data for an instance of 'alg'.
%
% Populate remaining fields of the toolbox data structure used by 'alg'
% function objects.
%
% DATA = ALG_INIT_DATA(DATA, X0, P0)
%
% DATA - Toolbox data structure.
% X0   - Initial solution guess for problem variables.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_init_data.m 2839 2015-03-05 17:09:01Z fschild $

xdim       = numel(x0);      % Number of problem variables
pdim       = numel(p0);      % Number of problem parameters
data.x_idx = (1:xdim)';      % Index set for problem variables
data.p_idx = xdim+(1:pdim)'; % Index set for problem parameters

end
