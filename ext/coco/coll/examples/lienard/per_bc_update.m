function data = per_bc_update(data, T, x0, x1, p) %#ok<INUSL>
%PER_BC_UPDATE   Update Poincare section.
%
% Use found trajectory end point at t=0 as new reference point and the
% corresponding vector field as new normal vector.
%
% see also per_bc and per_bc_DFDX

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: per_bc_update.m 2839 2015-03-05 17:09:01Z fschild $

n = numel(x0);
q = numel(p);

data.x0 = x0;
data.f0 = data.fhan(x0,p)';
data.J  = [sparse(n,1), speye(n,n), -speye(n,n), sparse(n,q);
           sparse(1,1), data.f0,    sparse(1,n), sparse(1,q)];

end
