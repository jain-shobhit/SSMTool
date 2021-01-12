function [fbc, Jbc] = torus_bc(data, T0, T, x0, x1, p)
%TORUS_BC   Torus boundary conditions.
%
% Trajectory end points lie on a curve on the invariant torus. The return
% map corresponds to identical times-of-flight and describes a rigid
% rotation.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: torus_bc.m 2839 2015-03-05 17:09:01Z fschild $

fbc = [T0; T-2*pi/p(1); data.F*x1-data.RF*x0; x0(2)];

nt = numel(T);
nx = numel(x0);
np = numel(p);

J1 = zeros(1,2*nt+2*nx+np);
J1(1,2*nt+2) = 1;

Jbc = [
  eye(nt), zeros(nt,nt+2*nx+np);
  zeros(nt), eye(nt), zeros(nt,2*nx), 2*pi/p(1)^2*ones(nt,1), zeros(nt,np-1);
  zeros(nx,2*nt), -data.RF, data.F, zeros(nx,np);
  J1];
end
