function prob = riess_start_1(prob, run, lab)
%RIESS_START_1   Append heteroclinic orbit problem constructed from stored data of periodic orbit.
%
% Construct an instance of 'po', append the corresponding variational
% zero problem, add segments in the stable manifold of the periodic orbit
% and the unstable manifold of the equilibrium at the origin, and introduce
% appropriate eigenspace and boundary conditions.
%
% PROB = RIESS_START_1(PROB, RUN, LAB)
%
% PROB - Continuation problem structure.
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: riess_start_1.m 2839 2015-03-05 17:09:01Z fschild $

[prob data sol] = povar_sol2orb(prob, '', run, lab); % Reconstruct orbit problem with variational zero problem

fdata = coco_get_func_data(prob, data.coll_id, 'data'); % Extract 'coll' data for periodic orbit
eps0  = [0.1; 0.1];
p0    = sol.x(fdata.p_idx); % Extract parameter values

s  = p0(1);
r  = p0(2);
t0 = 0; % Zero-length initial segment
x0 = eps0(1)*[(1-s+sqrt((1-s)^2+4*r*s))/2/r; 1; 0]'; % Point along unstable eigenvector at equilibrium
coll_args = {fdata.fhan, fdata.dfdxhan, fdata.dfdphan, t0, x0, p0};
prob = coll_isol2seg(prob, 'col1', coll_args{:}); % Append 'coll' instance

M      = reshape(sol.x(data.ubp_idx), data.u_shp); % Fundamental solution
M1     = M(data.M1_idx,:);                         % Jacobian of time-T flow
[v, d] = eig(M1);                                  % Floquet multipliers and eigenspaces
ind    = find(abs(diag(d))<1);
vec0   = -v(:,ind);                                % Stable eigenvector
lam0   = d(ind,ind);                               % Stable eigenvalue

t0 = 0; % Zero-length initial segment
x0 = sol.x(fdata.xbp_idx(end-data.dim+1:end))'+eps0(2)*vec0'; % Point along stable manifold of last point on orbit
coll_args = {fdata.fhan, fdata.dfdxhan, fdata.dfdphan, t0, x0, p0};
prob = coll_isol2seg(prob, 'col2', coll_args{:}); % Append 'coll' instance

data.nrm = [0, -1, 1]/sqrt(2); % Normal to hyperplane
data.pt0 = [20; 20; 30];       % Base point for hyperplane

prob = riess_close_het_1(prob, data, vec0, lam0, eps0); % Append eigenspace and boundary conditions

end
