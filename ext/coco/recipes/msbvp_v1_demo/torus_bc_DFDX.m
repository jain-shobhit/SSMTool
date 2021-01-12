function Jbc = torus_bc_DFDX(data, T, x0, x1, p)
%TORUS_BC_DFDX   Linearization of torus boundary conditions.
%
% Boundary conditions are linear and linearization is constant.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: torus_bc_DFDX.m 2839 2015-03-05 17:09:01Z fschild $

nt = numel(T);
nx = numel(x0);
np = numel(p);

J1 = zeros(2,nt+2*nx+np);
J1(1,nt+2) = 1;
J1(2,nt+[1 4]) = [-1 1];

Jbc = [eye(nt), zeros(nt,2*nx+np-1), -ones(nt,1);
     zeros(nx,nt), -data.RF, data.F, zeros(nx,np);
     J1];

end
