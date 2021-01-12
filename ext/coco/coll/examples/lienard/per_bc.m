function [fbc, Jbc] = per_bc(data, T, x0, x1, p) %#ok<INUSD,INUSL>
%PER_BC   Periodic boundary conditions on moving Poincare section and linearization
%
% Trajectory end point lies on hyperplane through reference point data.x0
% and with normal vector data.f0.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: per_bc.m 2839 2015-03-05 17:09:01Z fschild $

fbc = [x0-x1; data.f0*(x0-data.x0)];
Jbc = data.J;

end
