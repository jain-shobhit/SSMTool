function fbc = torus_bc(data, T, x0, x1, p)
%TORUS_BC   Torus boundary conditions.
%
% Trajectory end points lie on a curve on the invariant torus. The return
% map corresponds to identical times-of-flight and describes a rigid
% rotation.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: torus_bc.m 2839 2015-03-05 17:09:01Z fschild $

  fbc = [T-p(4); data.F*x1-data.RF*x0; x0(2); x0(4)-x0(1)];
end
