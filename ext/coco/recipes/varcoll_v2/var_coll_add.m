function prob = var_coll_add(prob, segoid, dfdxdx, dfdxdp)
%VAR_COLL_ADD   Append an instance of 'varcoll' to problem.
%
% Embed the variational problem associated with an instance of 'coll' as a
% zero problem in an existing continuation problem structure.
%
% PROB = VAR_COLL_ADD(PROB, SEGOID, DFDXDX, DFDXDP)
%
% PROB   - Continuation problem structure.
% SEGOID - 'coll' object instance identifier.
% DFDXDX - Function handle for vectorized encoding of 3d array of first
%          partials of vector-field Jacobian with respect to problem variables.
% DFDXDP - Function handle for vectorized encoding of 3d array of first
%          partials of vector-field Jacobian with respect to problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: var_coll_add.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(segoid, 'var');
data.coll_id   = coco_get_id(segoid, 'coll'); % 'coll' object instance identifier
data.dfdxdxhan = dfdxdx; % Second partials
data.dfdxdphan = dfdxdp; % Second partials

data = var_coll_init_data(prob, data); % Build toolbox data
M0   = var_coll_init_sol(prob, data);  % Build initial solution guess
uidx = coco_get_func_data(prob, data.coll_id, 'uidx'); % Extract 'coll' context-dependent index array
prob = coco_add_func(prob, tbid, @var_coll_F, @var_coll_DFDU, ...
  data, 'zero', 'uidx', uidx, 'u0', M0); % Append variational zero problem
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data y] = var_coll_F(prob, data, u)
%VAR_COLL_F   Variational collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.

fdata = coco_get_func_data(prob, data.coll_id, 'data'); % Extract 'coll' toolbox data

ubp = u(data.ubp_idx);
T   = u(fdata.T_idx);   % Extract interval length
x   = u(fdata.xbp_idx); % Extract basepoint values
p   = u(fdata.p_idx);   % Extract problem parameters

Mbp = reshape(ubp, data.u_shp);        % Fundamental solution at basepoints
xx  = reshape(fdata.W*x, fdata.x_shp); % Trajectory values at collocation nodes
pp  = repmat(p, fdata.p_rep);

ode = fdata.dfdxhan(xx, pp);
ode = sparse(fdata.dxrows, fdata.dxcols, ode(:));
ode = (0.5*T/fdata.coll.NTST)*ode*fdata.W-fdata.Wp;
ode = ode*Mbp;            % Collocation conditions
cnt = fdata.Q*Mbp;        % Continuity conditions
bcd = data.R*Mbp-data.Id; % Initial conditions

y = [ode(:); cnt(:); bcd(:)];

end

function [data J] = var_coll_DFDU(prob, data, u)
%VAR_COLL_DFDU   Linearization of variational collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.

fdata = coco_get_func_data(prob, data.coll_id, 'data'); % Extract 'coll' toolbox data

NTST = fdata.coll.NTST; % Number of mesh intervals

ubp = u(data.ubp_idx);
T   = u(fdata.T_idx);   % Extract interval length
x   = u(fdata.xbp_idx); % Extract basepoint values
p   = u(fdata.p_idx);   % Extract problem parameters

Mbp = reshape(ubp, data.u_shp);        % Fundamental solution at basepoints
xx  = reshape(fdata.W*x, fdata.x_shp); % Trajectory values at collocation nodes
pp  = repmat(p, fdata.p_rep);

dfdx = fdata.dfdxhan(xx, pp);
dfdx = sparse(fdata.dxrows, fdata.dxcols, dfdx(:));

dfdxdx = data.dfdxdxhan(xx, pp);
dfdxdx = sparse(data.dxdxrows1, data.dxdxcols1, dfdxdx(:));
dxode  = dfdxdx*fdata.W;
dxode  = sparse(data.dxdxrows2, data.dxdxcols2, dxode(:))*fdata.W*Mbp;
dxode  = (0.5*T/NTST)*sparse(data.dxdxrows3, data.dxdxcols3, dxode(:)); % W.r.t. xbp

dTode  = (0.5/NTST)*dfdx*fdata.W*Mbp; % W.r.t. T

dfdxdp = data.dfdxdphan(xx, pp);
dpode  = sparse(data.dxdprows, data.dxdpcols, dfdxdp(:));
dpode  = dpode*kron(speye(fdata.pdim), fdata.W*Mbp);
dpode  = (0.5*T/NTST)*reshape(dpode, data.dxdp_shp); % W.r.t. p

duode  = kron(data.Id, (0.5*T/NTST)*dfdx*fdata.W-fdata.Wp); % W.r.t. Mbp

J = [dxode, dTode(:), dpode, duode; data.jac];

end
