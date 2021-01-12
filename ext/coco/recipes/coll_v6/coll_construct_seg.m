function prob = coll_construct_seg(prob, tbid, data, sol)
%COLL_CONSTRUCT_SEG   Append an instance of 'coll' to problem.
%
% Add collocation and continuity conditions, monitor functions that
% evaluate to the problem parameters, and corresponding inactive
% continuation parameters. Append monitor function for discretization error
% estimate and terminate when this exceeds a given tolerance.
%
% Differs from coll_v5 by including a nonembedded continuation parameter
% for the number of discretization intervals.
%
% PROB = COLL_CONSTRUCT_SEG(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_seg.m 2923 2015-10-27 20:53:08Z hdankowicz $

data.tbid = tbid;
data = coco_func_data(data); % Convert to func_data class for shared access
prob = coco_add_func(prob, tbid, @coll_F, @coll_DFDU, data, 'zero', ...
  'u0', sol.u, 'remesh', @coll_remesh);
uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  prob = coco_add_pars(prob, fid, uidx(data.maps.p_idx), data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');
efid = coco_get_id(tbid, {'err' 'err_TF' 'NTST'});
prob = coco_add_chart_data(prob, tbid, struct(), struct()); % Allocate chart data
prob = coco_add_func(prob, efid{1}, @coll_err, data, ...
  'regular', efid, 'uidx', uidx(data.maps.xbp_idx), ...
  'remesh', @coll_err_remesh, 'passChart'); % Monitor discretization error estimate with remesh support and chart data
prob = coco_add_event(prob, 'MXCL', 'MX', efid{2}, '>', 1); % MXCL - terminal event type

end

function [data y] = coll_F(prob, data, u)
%COLL_F   Collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.
%
% Identical to coll_v5.

maps = data.maps;
mesh = data.mesh;

x = u(maps.xbp_idx); % Extract basepoint values
T = u(maps.T_idx);   % Extract interval length
p = u(maps.p_idx);   % Extract problem parameters

xx = reshape(maps.W*x, maps.x_shp); % Values at collocation nodes
pp = repmat(p, maps.p_rep);

ode = mesh.fka.*data.fhan(xx, pp);
ode = (0.5*T/maps.NTST)*ode(:)-maps.Wp*x; % Collocation conditions
cnt = maps.Q*x;                           % Continuity conditions

y = [ode; cnt];

end

function [data J] = coll_DFDU(prob, data, u)
%COLL_DFDU   Linearization of collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.
%
% Identical to coll_v5.

maps = data.maps;
mesh = data.mesh;

x = u(maps.xbp_idx); % Extract basepoint values
T = u(maps.T_idx);   % Extract interval length
p = u(maps.p_idx);   % Extract problem parameters

xx = reshape(maps.W*x, maps.x_shp); % Values at collocation nodes
pp = repmat(p, maps.p_rep);

if isempty(data.dfdxhan)
  dxode = coco_ezDFDX('f(x,p)v', data.fhan, xx, pp);
else
  dxode = data.dfdxhan(xx, pp);
end
dxode = mesh.dxka.*dxode;
dxode = sparse(maps.dxrows, maps.dxcols, dxode(:));
dxode = (0.5*T/maps.NTST)*dxode*maps.W-maps.Wp; % W.r.t. basepoint values

dTode = mesh.fka.*data.fhan(xx, pp);
dTode = (0.5/maps.NTST)*dTode(:); % W.r.t. interval length

if isempty(data.dfdphan)
  dpode = coco_ezDFDP('f(x,p)v', data.fhan, xx, pp);
else
  dpode = data.dfdphan(xx, pp);
end
dpode = mesh.dpka.*dpode;
dpode = sparse(maps.dprows, maps.dpcols, dpode(:));
dpode = (0.5*T/maps.NTST)*dpode; % W.r.t. interval length

J = [dxode dTode dpode; maps.Q sparse(maps.Qnum,1+maps.pdim)];

end

function [prob stat xtr] = coll_remesh(prob, data, chart, ub, Vb)
%COLLREMESH   Remesh function for equidistribution of error estimate.
%
% Update warping coefficients to reflect an equidistribution of the
% interpolation error estimate.
%
% Differs from coll_v5 by providing support for varying discretization
% order.

int  = data.int;
maps = data.maps;
mesh = data.mesh;
coll = data.coll;

xtr  = maps.xtr; % Current invariant elements

u = ub(maps.xbp_idx);   % Extract current basepoint values
V = Vb(maps.xbp_idx,:); % Extract current tangent vector and tangent space

cp   = reshape(maps.Wm*u, [int.dim maps.NTST]);
cp   = sqrt(sum(cp.^2,1));
cp   = nthroot(cp, int.NCOL); % Piecewise constant value of dF on each interval
cpmn = nthroot(0.5*(coll.TOLINC+coll.TOLDEC)/maps.wn, int.NCOL);

ka = mesh.ka;
s  = data.coll.SAD; % Equidistributed when s=1, uniform when s=0
F  = [0 cumsum((1-s)*cpmn*ka + s*cp, 2)]; % Cumulative error estimate
t  = [0 cumsum(ka, 2)];

NTSTi = min(ceil(maps.NTST*1.1+1),   coll.NTSTMX);
NTSTd = max(ceil(maps.NTST/1.025-1), coll.NTSTMN);
cdata = coco_get_chart_data(chart, data.tbid);
err   = cdata.err(1);
maps2 = maps;
if err>coll.TOLINC && maps.NTST~=NTSTi % If above upper bound on adaptation window
  data.coll.NTST = NTSTi;
  maps2 = coll_maps(int, NTSTi, maps.pdim); % Update toolbox data
elseif err<coll.TOLDEC && maps.NTST~=NTSTd % If below lower bound on adaptation window
  data.coll.NTST = NTSTd;
  maps2 = coll_maps(int, NTSTd, maps.pdim); % Update toolbox data
end
xtr(maps.xtrend) = maps2.xtr(maps2.xtrend); % Translation table
th = linspace(0, F(end), maps2.NTST+1);
tt = interp1(F, t*(maps2.NTST/t(end)), th, 'pchip'); % Distribute error estimate across new mesh

mesh2 = coll_mesh(int, maps2, tt); % Update mesh-dependent data

tbp = mesh.tbp(maps.tbp_idx); % New basepoints
xbp = reshape(u, maps.xbp_shp);
xbp = xbp(:, maps.tbp_idx);
x0  = interp1(tbp', xbp', mesh2.tbp, 'pchip')'; % Update basepoint values
u1  = x0(:);

V1 = zeros(numel(mesh2.tbp)*int.dim, size(V,2));
for i=1:size(V,2)
  vbp = reshape(V(:,i), maps.xbp_shp);
  vbp = vbp(:, maps.tbp_idx);
  v0  = interp1(tbp', vbp', mesh2.tbp, 'pchip')'; % Update tangent vectors and tangent space
  V1(:,i) = v0(:);
end

data.maps = maps2;
data.mesh = mesh2;

ua = [u1; ub(maps.Tp_idx)];
Va = [V1; Vb(maps.Tp_idx,:)];

prob = coco_change_func(prob, data, 'u0', ua, 'vecs', Va); % Update data and solution

stat = 'success';

end

function [data chart y] = coll_err(prob, data, chart, u)
%COLL_ERR   Evaluate estimate of approximation eror.
%
% Estimate discretization error using the highest-order coefficients of the
% interpolating Lagrange polynomials. Return estimated error and scaled by
% error tolerance.
%
% Differs from coll_v5 by including the number of discretization intervals
% in the y return argument.

cdata = coco_get_chart_data(chart, data.tbid);
if isfield(cdata, 'err')
  y = cdata.err;
else
  int  = data.int;
  maps = data.maps;
  
  cp = reshape(maps.Wm*u, [int.dim maps.NTST]);
  y  = maps.wn*max(sqrt(sum(cp.^2,1)));
  y  = [y; y/data.coll.TOL; data.maps.NTST];
  cdata.err = y;
  chart = coco_set_chart_data(chart, data.tbid, cdata);
end

end

function [prob stat xtr] = coll_err_remesh(prob, data, chart, ub, Vb)
%COLL_ERR_REMESH   Update dependency index set after remeshing.
%
% Associate a remesh action with the coll_err function.
%
% Identical to coll_v5.

maps = data.maps;

xtr    = []; % No invariant indices
uidx   = coco_get_func_data(prob, data.tbid, 'uidx');
prob   = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx)); % Update data and solution
stat = 'success';

end
