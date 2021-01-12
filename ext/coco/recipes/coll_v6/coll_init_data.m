function data = coll_init_data(data, t0, x0, p0)
%COLL_INIT_DATA   Initialize toolbox data for an instance of 'coll'.
%
% Populate remaining fields of the toolbox data structure used by 'coll'
% function objects.
%
% Identical to coll_v5.
%
% DATA = COLL_INIT_DATA(DATA, T0, X0, P0)
%
% DATA - Toolbox data structure.
% T0   - Initial temporal mesh.
% X0   - Initial solution guess for discretized trajectory.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_init_data.m 2839 2015-03-05 17:09:01Z fschild $

NCOL = data.coll.NCOL; % Degree of polynomial interpolants
dim  = size(x0,2);     % State-space dimension
data.int = coll_interval(NCOL, dim); % Props independent of NTST and mesh distribution

NTST = data.coll.NTST; % Number of mesh intervals
pdim = numel(p0);      % Number of problem parameters
data.maps = coll_maps(data.int, NTST, pdim); % Props depend on NTST, but not on mesh distribution

t  = linspace(0, NTST, numel(t0));
tt = interp1(t, t0, 0:NTST, 'linear');
tt = tt*(NTST/tt(end)); % Initial temporal mesh
data.mesh = coll_mesh(data.int, data.maps, tt); % Props depend on mesh distribution

end
