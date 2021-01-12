function Jbc = per_bc_DFDX(data, T, x0, x1, p)
%PER_BC_DFDX   Linearization of periodic boundary conditions on moving Poincare section
%
% Boundary conditions are linear and linearization is constant.
%
% see also per_bc_update

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: per_bc_DFDX.m 2839 2015-03-05 17:09:01Z fschild $

  Jbc = data.J;
end
