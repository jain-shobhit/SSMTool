function data = alg_init_data(data, x0, p0)
%ALG_INIT_DATA   Initialize toolbox data for an instance of 'alg'.
%
% Populate remaining fields of the toolbox data structure used by 'alg'
% function objects.
%
% Differs from alg_v4 by the construction of toolbox data particular to the
% bordered matrix algorithm used for fold detection.
%
% DATA = ALG_INIT_DATA(DATA, X0, P0)
%
% DATA - Toolbox data structure.
% X0   - Initial solution guess for problem variables.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_init_data.m 2955 2017-01-10 15:19:26Z hdankowicz $

xdim       = numel(x0);      % Number of problem variables
pdim       = numel(p0);      % Number of problem parameters
data.x_idx = (1:xdim)';      % Index set for problem variables
data.p_idx = xdim+(1:pdim)'; % Index set for problem parameters

Jx         = alg_fhan_DFDX(data, x0, p0);
[data.b, ~, data.c] = svds(Jx,1,'smallest'); % Bordering vectors
data.rhs   = [zeros(xdim,1); 1];    % Right-hand side

end
