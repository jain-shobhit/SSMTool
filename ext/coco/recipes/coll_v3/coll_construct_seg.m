function prob = coll_construct_seg(prob, tbid, data, sol)
%COLL_CONSTRUCT_SEG   Append an instance of 'coll' to problem.
%
% Add collocation and continuity conditions, monitor functions that
% evaluate to the problem parameters, and corresponding inactive
% continuation parameters. Append monitor function for discretization error
% estimate and terminate when this exceeds a given tolerance.
%
% Differs from coll_v1 by providing initial support for adaptive changes to
% discretization order.
%
% PROB = COLL_CONSTRUCT_SEG(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_seg.m 2839 2015-03-05 17:09:01Z fschild $

prob = coco_add_func(prob, tbid, @coll_F, @coll_DFDU, data, 'zero', ...
  'u0', sol.u);
uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  prob = coco_add_pars(prob, fid, uidx(data.p_idx), data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');
efid = coco_get_id(tbid, {'err' 'err_TF'});
prob = coco_add_func(prob, efid{1}, @coll_err, data, ...
  'regular', efid, 'uidx', uidx(data.xbp_idx)); % Monitor discretization error estimate
prob = coco_add_event(prob, 'MXCL', 'MX', efid{2}, '>', 1); % MXCL - terminal event type

end

function [data y] = coll_F(prob, data, u)
%COLL_F   Collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.
%
% Identical to coll_v1.

x = u(data.xbp_idx); % Extract basepoint values
T = u(data.T_idx);   % Extract interval length
p = u(data.p_idx);   % Extract problem parameters

xx = reshape(data.W*x, data.x_shp); % Values at collocation nodes
pp = repmat(p, data.p_rep);

ode = data.fhan(xx, pp);
ode = (0.5*T/data.coll.NTST)*ode(:)-data.Wp*x; % Collocation conditions
cnt = data.Q*x;                                % Continuity conditions

y = [ode; cnt];

end

function [data J] = coll_DFDU(prob, data, u)
%COLL_DFDU   Linearization of collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.
%
% Minor difference from coll_v1 in use of Qnum instead of dTpcnt.

x = u(data.xbp_idx); % Extract basepoint values
T = u(data.T_idx);   % Extract interval length
p = u(data.p_idx);   % Extract problem parameters

xx = reshape(data.W*x, data.x_shp); % Values at collocation nodes
pp = repmat(p, data.p_rep);

if isempty(data.dfdxhan)
  dxode = coco_ezDFDX('f(x,p)v', data.fhan, xx, pp);
else
  dxode = data.dfdxhan(xx, pp);
end
dxode = sparse(data.dxrows, data.dxcols, dxode(:));
dxode = (0.5*T/data.coll.NTST)*dxode*data.W-data.Wp; % W.r.t. basepoint values

dTode = data.fhan(xx, pp);
dTode = (0.5/data.coll.NTST)*dTode(:); % W.r.t. interval length

if isempty(data.dfdphan)
  dpode = coco_ezDFDP('f(x,p)v', data.fhan, xx, pp);
else
  dpode = data.dfdphan(xx, pp);
end
dpode = sparse(data.dprows, data.dpcols, dpode(:));
dpode = (0.5*T/data.coll.NTST)*dpode; % W.r.t. problem parameters

J = [dxode dTode dpode; data.Q sparse(data.Qnum,1+data.pdim)];

end

function [data y] = coll_err(prob, data, u)
%COLL_ERR   Evaluate estimate of approximation eror.
%
% Estimate discretization error using the highest-order coefficients of the
% interpolating Lagrange polynomials. Return estimated error and scaled by
% error tolerance.

cp = reshape(data.Wm*u, [data.dim data.coll.NTST]);
y  = data.wn*max(sqrt(sum(cp.^2,1)));
y  = [y; y/data.coll.TOL];

end
