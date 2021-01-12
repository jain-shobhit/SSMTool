function prob = coll_construct_seg(prob, tbid, data, sol)
%COLL_CONSTRUCT_SEG   Append an instance of 'coll' to problem.
%
% Add collocation and continuity conditions, monitor functions that
% evaluate to the problem parameters, and corresponding inactive
% continuation parameters. Append monitor function for discretization error
% estimate and terminate when this exceeds a given tolerance.
%
% Identical to coll_v4.
%
% PROB = COLL_CONSTRUCT_SEG(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_seg.m 2839 2015-03-05 17:09:01Z fschild $

data.tbid = tbid;
data = coco_func_data(data); % Convert to func_data class for shared access
prob = coco_add_func(prob, tbid, @coll_F, @coll_DFDU, data, 'zero', ...
  'u0', sol.u);
uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  prob = coco_add_pars(prob, fid, uidx(data.maps.p_idx), data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');
efid = coco_get_id(tbid, {'err' 'err_TF'});
prob = coco_add_func(prob, efid{1}, @coll_err, data, ...
  'regular', efid, 'uidx', uidx(data.maps.xbp_idx)); % Monitor discretization error estimate
prob = coco_add_event(prob, 'MXCL', 'MX', efid{2}, '>', 1); % MXCL - terminal event type
prob = coco_add_slot(prob, tbid, @coll_update, data, 'update'); % Update reference warping coefficients
prob = coco_add_slot(prob, tbid, @coll_update_h, data, 'update_h'); % Update mesh condition step size

end

function [data y] = coll_F(prob, data, u)
%COLL_F   Collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.
%
% Identical to coll_v4.

maps = data.maps;

x  = u(maps.xbp_idx); % Extract basepoint values
T  = u(maps.T_idx);   % Extract interval length
p  = u(maps.p_idx);   % Extract problem parameters
ka = u(maps.ka_idx);  % Extract warping coefficients
la = u(maps.la_idx);  % Extract Lagrange multiplier

fka = ka(maps.fka_idx);
xx  = reshape(maps.W*x, maps.x_shp); % Values at collocation nodes
pp  = repmat(p, maps.p_rep);

ode = fka.*data.fhan(xx, pp);
ode = (0.5*T/maps.NTST)*ode(:)-maps.Wp*x; % Collocation conditions
cnt = maps.Q*x;                           % Continuity conditions

v   = reshape(x, maps.v_shp);
msh = v(:,1:end-1,:)-v(:,2:end,:);
msh = squeeze(sum(sqrt(sum(msh.*msh,1)),2));
msh = [ka-la*data.ka+data.h*(msh-data.mean_msh); sum(ka)-maps.NTST]; % Mesh conditions and conservation condition

y = [ode; cnt; msh];

end

function [data J] = coll_DFDU(prob, data, u)
%COLL_DFDU   Linearization of collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.
%
% Identical to coll_v4.

maps = data.maps;
int  = data.int;

x  = u(maps.xbp_idx); % Extract basepoint values
T  = u(maps.T_idx);   % Extract interval length
p  = u(maps.p_idx);   % Extract problem parameters
ka = u(maps.ka_idx);  % Extract warping coefficients

fka  = ka(maps.fka_idx);
dxka = ka(maps.dxka_idx);
dpka = ka(maps.dpka_idx);

xx = reshape(maps.W*x, maps.x_shp); % Values at collocation nodes
pp = repmat(p, maps.p_rep);

if isempty(data.dfdxhan)
  dxode = coco_ezDFDX('f(x,p)v', data.fhan, xx, pp);
else
  dxode = data.dfdxhan(xx, pp);
end
dxode = dxka.*dxode;
dxode = sparse(maps.dxrows, maps.dxcols, dxode(:));
dxode = (0.5*T/maps.NTST)*dxode*maps.W-maps.Wp; % W.r.t. basepoint values

ode   = data.fhan(xx, pp);
dTode  = fka.*ode;
dTode  = (0.5/maps.NTST)*dTode(:); % W.r.t. interval length
dkaode = (0.5*T/maps.NTST)*ode;
dkaode = sparse(maps.karows, maps.fka_idx, dkaode(:)); % W.r.t. warping coefficients

if isempty(data.dfdphan)
  dpode = coco_ezDFDP('f(x,p)v', data.fhan, xx, pp);
else
  dpode = data.dfdphan(xx, pp);
end
dpode = dpka.*dpode;
dpode = sparse(maps.dprows, maps.dpcols, dpode(:));
dpode = (0.5*T/maps.NTST)*dpode; % W.r.t. problem parameters

J1 = [dxode dTode dpode dkaode; ...
  maps.Q sparse(maps.Qnum, 1+maps.pdim+maps.NTST)];

v    = reshape(x, maps.v_shp);
du   = v(:,1:end-1,:)-v(:,2:end,:);
dsq  = 1./sqrt(sum(du.*du,1));
dsq  = du.*repmat(dsq, [int.dim 1 1]);

df1  = dsq(:,1,:);
df2  = -dsq(:,1:end-1,:)+dsq(:,2:end,:);
df3  = -dsq(:,end,:);
df   = data.h*cat(2,df1,df2,df3);

rows = reshape(1:maps.NTST, [1 1 maps.NTST]);
rows = repmat(rows, [int.dim int.NCOL+1 1]);
cols = reshape(maps.xbp_idx, maps.v_shp);
vals = df;

rows = [rows(:); (1:maps.NTST)'];
cols = [cols(:); maps.ka_idx];
vals = [vals(:); ones(maps.NTST,1)];

rows = [rows(:); (1:maps.NTST)'];
cols = [cols(:); maps.la_idx*ones(maps.NTST,1)];
vals = [vals(:); -data.ka];

rows = [rows(:); (maps.NTST+1)*ones(maps.NTST,1)];
cols = [cols(:); maps.ka_idx];
vals = [vals(:); ones(maps.NTST,1)];

J3 = sparse(rows, cols, vals); % Of mesh and conservation conditions

J = sparse([J1 zeros(int.dim*(int.NCOL+1)*maps.NTST-int.dim,1); J3]);

end

function [data y] = coll_err(prob, data, u)
%COLL_ERR   Evaluate estimate of approximation eror.
%
% Estimate discretization error using the highest-order coefficients of the
% interpolating Lagrange polynomials. Return estimated error and scaled by
% error tolerance.
%
% Identical to coll_v4.

int  = data.int;
maps = data.maps;
  
cp = reshape(maps.Wm*u, [int.dim maps.NTST]);
y  = maps.wn*max(sqrt(sum(cp.^2,1)));
y  = [y; y/data.coll.TOL];

end

function data = coll_update(prob, data, cseg, varargin)
%COLL_UPDATE   Update reference warping coefficients and arclengths.
%
% Slot function for updating reference values appearing in the
% comoving-mesh conditions.
%
% Differs from coll_v4 by inclusion of scaling factor for step size.

maps = data.maps;
int  = data.int;

data.h = data.coll.hfac*cseg.prcond.h; % Scaled projection condition step size

base_chart = cseg.src_chart;
uidx       = coco_get_func_data(prob, data.tbid, 'uidx');
u          = base_chart.x(uidx);
data.ka    = u(maps.ka_idx);
tmi        = [0 cumsum(data.ka')];
data.mesh  = coll_mesh(data.int, data.maps, tmi); % Reference warp

v   = reshape(u(maps.xbp_idx), [int.dim int.NCOL+1 maps.NTST]);
msh = v(:,1:end-1,:)-v(:,2:end,:);
msh = squeeze(sum(sqrt(sum(msh.*msh,1)),2));
data.mean_msh = mean(msh); % Average interval arclength

end

function data = coll_update_h(prob, data, h, varargin)
%COLL_UPDATE_H   Update step size of mesh condition.
% 
% Slot function: Use projection condition step size to update step size of
% comoving-mesh condition.
%
% Differs from coll_v4 by inclusion of scaling factor for step size.

data.h = data.coll.hfac*h;

end
