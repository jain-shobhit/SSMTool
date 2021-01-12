function [prob, data] = coll_construct_seg(prob, data, sol, CacheJacobian)
%COLL_CONSTRUCT_SEG   Add COLL zero problem.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_seg.m 3006 2017-04-26 12:14:37Z hdankowicz $

[data, sol] = init_data(data, sol);
seg  = data.coll_seg;
fid  = seg.fid;

if CacheJacobian
  data.sh.coll_Jx = [];
  prob = coco_add_func(prob, fid, @FDF1, data, 'zero', ...
    'u0', sol.u0, 't0', sol.t0, 'remesh', @remesh, 'F+DF');%, ...
    %'prop', 'weights', {seg.maps.xbp_idx, (0.5/seg.maps.NTST)*seg.maps.W'*seg.mesh.kas2*seg.mesh.wts2*seg.maps.W});
  data = data.protect('coll_Jx');
  data.no_save = [ data.no_save { 'coll_Jx' } ];
else
  prob = coco_add_func(prob, fid, @FDF2, data, 'zero', ...
    'u0', sol.u0, 't0', sol.t0, 'remesh', @remesh, 'F+DF');%, ...
    %'prop', 'weights', {seg.maps.xbp_idx, (0.5/seg.maps.NTST)*seg.maps.W'*seg.mesh.kas2*seg.mesh.wts2*seg.maps.W});
end

if data.ode.autonomous
  uidx = coco_get_func_data(prob, fid, 'uidx');
  tfid = coco_get_id(fid, 'T0');
  prob = coco_add_pars(prob, tfid, uidx(seg.maps.T0_idx), tfid);
end

if ~isempty(data.pnames)
  uidx = coco_get_func_data(prob, fid, 'uidx');
  pfid = coco_get_id(fid, 'pars');
  prob = coco_add_pars(prob, pfid, uidx(seg.maps.p_idx), data.pnames);
end

end

function [data, sol] = init_data(data, sol)

t0 = sol.tbp-sol.T0;

NCOL = data.coll.NCOL;
seg.int  = coll_interval(NCOL, data.xdim);

NTST = data.coll.NTST;
seg.maps = coll_maps(seg.int, NTST, data.pdim);

if abs(sol.T)>eps
  t  = linspace(0, NTST, numel(t0));
  tt = interp1(t, t0, 0:NTST, 'linear');
  tt = tt*(NTST/tt(end));
else
  tt = 0:NTST;
end
seg.mesh = coll_mesh(seg.int, seg.maps, tt);

seg.fid  = coco_get_id(data.oid, 'coll');
data.coll_seg = seg;

x0 = sol.xbp;

if abs(sol.T)>eps
  t0 = t0/sol.T;
  x0 = interp1(t0, x0, seg.mesh.tbp)';
else
  x0 = repmat(x0(1,:), size(seg.mesh.tbp))';
end

sol.u0 = [x0(:); sol.T0; sol.T; sol.p];
  
if ~isempty(sol.t0)
  x0_t = sol.xbp_t0;
  if abs(sol.T)>eps
    x0_t = interp1(t0, x0_t, seg.mesh.tbp)';
  else
    x0_t = repmat(x0_t(1,:), size(seg.mesh.tbp))';
  end
  sol.t0 = [x0_t(:); sol.T0_t0; sol.T_t0; sol.p_t0];
end

end

function [data, y, J] = FDF1(~, data, u) 
%FDF1   COLL zero problem with cacheing of derivatives.

pr = data.pr;
seg  = pr.coll_seg;
maps = seg.maps;
mesh = seg.mesh;

x  = u(maps.xbp_idx);
T0 = u(maps.T0_idx);
T  = u(maps.T_idx);
p  = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
pcn = repmat(p, maps.p_rep);
fcn = pr.ode_F(pr, T0+T*mesh.tcn', xcn, pcn);
ode = mesh.fka.*fcn;
ode = (0.5*T/maps.NTST)*ode(:)-maps.Wp*x;
cnt = maps.Q*x;

y = [ode; cnt];

fdxcn = pr.ode_DFDX(pr, T0+T*mesh.tcn', xcn, pcn);
dxode = mesh.fdxka.*fdxcn;
dxode = sparse(maps.fdxrows, maps.fdxcols, dxode(:));
dxode = (0.5*T/maps.NTST)*dxode*maps.W-maps.Wp;

data.coll_Jx = dxode;

if nargout>=3
  fdtcn = pr.ode_DFDT(pr, T0+T*mesh.tcn', xcn, pcn);
  dT0ode = mesh.fka.*(T*fdtcn);
  dT0ode = (0.5/maps.NTST)*dT0ode(:);
  dTode = mesh.fka.*(fcn+T*fdtcn.*repmat(mesh.tcn',[pr.xdim 1]));
  dTode = (0.5/maps.NTST)*dTode(:);
 
  fdpcn = pr.ode_DFDP(pr, T0+T*mesh.tcn', xcn, pcn);
  dpode = mesh.fdpka.*fdpcn;
  dpode = sparse(maps.fdprows, maps.fdpcols, dpode(:));
  dpode = (0.5*T/maps.NTST)*dpode;
  
  J = [dxode dT0ode dTode dpode; maps.Q sparse(maps.Qnum,2+maps.pdim)];
end

end

function [data, y, J] = FDF2(prob, data, u) 
%FDF   COLL zero problem without cacheing of derivatives.

pr = data.pr;
seg  = pr.coll_seg;
maps = seg.maps;
mesh = seg.mesh;

x  = u(maps.xbp_idx);
T0 = u(maps.T0_idx);
T  = u(maps.T_idx);
p  = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
pcn = repmat(p, maps.p_rep);
fcn = pr.ode_F(pr, T0+T*mesh.tcn', xcn, pcn);
ode = mesh.fka.*fcn;
ode = (0.5*T/maps.NTST)*ode(:)-maps.Wp*x;
cnt = maps.Q*x;

y = [ode; cnt];

if nargout>=3
  fdxcn = pr.ode_DFDX(pr, T0+T*mesh.tcn', xcn, pcn);
  dxode = mesh.fdxka.*fdxcn;
  dxode = sparse(maps.fdxrows, maps.fdxcols, dxode(:));
  dxode = (0.5*T/maps.NTST)*dxode*maps.W-maps.Wp;
  
  fdtcn = pr.ode_DFDT(pr, T0+T*mesh.tcn', xcn, pcn);
  dT0ode = mesh.fka.*(T*fdtcn);
  dT0ode = (0.5/maps.NTST)*dT0ode(:);
  dTode = mesh.fka.*(fcn+T*fdtcn.*repmat(mesh.tcn',[pr.xdim 1]));
  dTode = (0.5/maps.NTST)*dTode(:);
  
  fdpcn = pr.ode_DFDP(pr, T0+T*mesh.tcn', xcn, pcn);
  dpode = mesh.fdpka.*fdpcn;
  dpode = sparse(maps.fdprows, maps.fdpcols, dpode(:));
  dpode = (0.5*T/maps.NTST)*dpode;
  
  J = [dxode dT0ode dTode dpode; maps.Q sparse(maps.Qnum,2+maps.pdim)];
end

end

function [prob, stat, xtr] = remesh(prob, data, chart, ub, Vb)
%COLL_REMESH   Remesh function for equidistribution of error estimate.
%
% Update warping coefficients to reflect an equidistribution of the
% interpolation error estimate.

pr = data.pr;
seg  = pr.coll_seg;
int  = seg.int;
maps = seg.maps;
mesh = seg.mesh;
coll = pr.coll;

xtr  = maps.xtr; % Current invariant elements

u = ub(maps.xbp_idx);   % Extract current basepoint values
V = Vb(maps.xbp_idx,:); % Extract current tangent vector and tangent space

cp   = reshape(maps.Wm*u, [int.dim maps.NTST]);
cp   = sqrt(sum(cp.^2,1));
cp   = nthroot(cp, int.NCOL); % Piecewise constant value of dF on each interval
cpmn = nthroot(0.5*(coll.TOLINC+coll.TOLDEC)/maps.wn, int.NCOL);

ka = mesh.ka;
s  = coll.SAD; % Equidistributed when s=1, uniform when s=0
F  = [0 cumsum((1-s)*cpmn*ka + s*cp, 2)]; % Cumulative error estimate
t  = [0 cumsum(ka, 2)];

NTSTi = min(ceil(maps.NTST*1.1+1),   coll.NTSTMX);
NTSTd = max(ceil(maps.NTST/1.025-1), coll.NTSTMN);
cdata = coco_get_chart_data(chart, seg.fid);
err   = cdata.err(1);
maps2 = maps;
if err>coll.TOLINC && maps.NTST~=NTSTi % If above upper bound on adaptation window
  pr.coll.NTST = NTSTi;
  maps2 = coll_maps(int, NTSTi, maps.pdim); % Update toolbox data
elseif err<coll.TOLDEC && maps.NTST~=NTSTd % If below lower bound on adaptation window
  pr.coll.NTST = NTSTd;
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

seg.maps = maps2;
seg.mesh = mesh2;
pr.coll_seg = seg;

ua = [u1; ub(maps.Tp_idx)];
Va = [V1; Vb(maps.Tp_idx,:)];

data.pr = pr;
prob = coco_change_func(prob, data, 'u0', ua, 'fdim', maps2.fdim, ...
  'vecs', Va);%, ...
    %'prop', 'weights', {seg.maps.xbp_idx, (0.5/seg.maps.NTST)*seg.maps.W'*seg.mesh.kas2*seg.mesh.wts2*seg.maps.W}); % Update data and solution

stat = 'success';

end
