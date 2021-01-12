function [prob, data] = coll_construct_var(prob, data)
%COLL_CONSTRUCT_VAR   Add variational problem as empty test function.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_var.m 2998 2017-01-29 16:25:42Z hdankowicz $

data = init_data(prob, data);
tst  = data.coll_tst;
fid  = tst.fid;

uidx = coco_get_func_data(prob, data.coll_seg.fid, 'uidx');
prob = coco_add_chart_data(prob, fid, [], []);
prob = coco_add_func(prob, fid, @var, data, 'regular', {}, ...
  'uidx', uidx, 'remesh', @remesh, 'passChart');
prob = coco_add_slot(prob, fid, @var_update, data, 'update');

end

function data = init_data(prob, data)

seg  = data.coll_seg;
dim  = data.xdim;
mesh = seg.mesh;
maps = seg.maps;

tst.fid  = coco_get_id(seg.fid, 'test');

tst.M0_idx = maps.xbp_idx(1:dim); % maps.x0_idx
tst.M1_idx = maps.xbp_idx(end-dim+(1:dim)); % maps.x1_idx
tst.row0   = [eye(dim) zeros(dim, maps.xbp_idx(end)-dim)];
tst.rhs0   = zeros(maps.xbp_idx(end)-dim, dim);

% tst.tbp     = mesh.tbp;
% tst.tbp_idx = maps.tbp_idx;
% tst.xbp_shp = maps.xbp_shp;

switch data.coll.method
  case '2I'
    wts     = kron(mesh.gwt, eye(dim));
    kas     = kron(mesh.gka, eye(dim));
    tst.B1  = (0.5/maps.NTST)*((wts.*kas)*maps.W); % Check kappas!
    tst.rhs = [2*eye(dim); zeros(maps.xbp_idx(end)-dim, dim)];
    tst.la1 = 0.5;
    tst.la2 = 0.5;
  case '3I'
    ze = sparse(dim, maps.xbp_idx(end)-2*dim);
    tst.B1  = [speye(dim) ze speye(dim)];
    tst.rhs = [3*eye(dim); zeros(maps.xbp_idx(end)-dim, dim)];
    tst.la1 = 2/3;
    tst.la2 = 1-tst.la1; % don't use 1/3 [roundoff errors]
  otherwise
    error('%s: unknown mode', mfilename);
end
tst.B2 = (0.5/data.coll.NTST)*maps.W'*mesh.wts2*mesh.kas2*maps.W; % Check kappas!

u = coco_get_func_data(prob, seg.fid, 'u0');

x  = u(maps.xbp_idx);
T0 = u(maps.T0_idx);
T  = u(maps.T_idx);
p  = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
pcn = repmat(p, maps.p_rep);

fdxcn = data.ode_DFDX(data, T0+T*mesh.tcn', xcn, pcn);
dxode = mesh.fdxka.*fdxcn;
dxode = sparse(maps.fdxrows, maps.fdxcols, dxode(:));
dxode = (0.5*T/maps.NTST)*dxode*maps.W-maps.Wp;

J0 = [dxode; maps.Q];

M0 = repmat(eye(dim), [data.coll.NTST*(data.coll.NCOL+1) 1]);
M  = [tst.B1+M0'*tst.B2; J0]\tst.rhs;
betas = linspace(0,1,max(1,data.coll.NBeta)+1);
for beta = betas(2:end)
  J1 = beta*J0;
  It = 1;
  while It<data.coll.NBItMX && norm((M(:)-M0(:)))>0.1
    M0 = tst.la1*M0 + tst.la2*M;
    J  = [tst.B1+M0'*tst.B2; J1];
    M  = J\tst.rhs;
    It = It+1;
  end
end
tst.M0  = tst.la1*M0 + tst.la2*M;
tst.row = tst.B1+M'*tst.B2;

data.coll_tst = tst;
data.no_save = [ data.no_save { 'coll_tst.B1' 'coll_tst.B2' ...
  'coll_tst.row' 'coll_tst.row0' 'coll_tst.M' 'coll_tst.M0' ...
  'coll_tst.rhs' 'coll_tst.rhs0' 'coll_tst.la1' 'coll_tst.la2' } ];

end

function [data, chart, y] = var(prob, data, chart, u) %#ok<INUSL>

pr = data.pr;
seg  = pr.coll_seg;
tst  = pr.coll_tst;
maps = seg.maps;
mesh = seg.mesh;

x    = u(maps.xbp_idx);
T0   = u(maps.T0_idx);
T    = u(maps.T_idx);
p    = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
pcn = repmat(p, maps.p_rep);

fdxcn = pr.ode_DFDX(pr, T0+T*mesh.tcn', xcn, pcn);
dfode = mesh.fdxka.*fdxcn;
dfode = sparse(maps.fdxrows, maps.fdxcols, dfode(:));
dfode = (0.5*T/maps.NTST)*dfode*maps.W-maps.Wp;

J     = [tst.row; dfode; maps.Q];
tst.M = J\tst.rhs;

data.coll_tst = tst;

cdata.M0 = tst.M(tst.M0_idx,:);
chart  = coco_set_chart_data(chart, tst.fid, cdata);

y = [];

end

function data = var_update(prob, data, cseg, varargin)

pr = data.pr;
seg  = pr.coll_seg;
tst  = pr.coll_tst;
maps = seg.maps;
mesh = seg.mesh;

uidx = coco_get_func_data(prob, seg.fid, 'uidx');
u    = cseg.src_chart.x(uidx);

x  = u(maps.xbp_idx);
T0 = u(maps.T0_idx); 
T  = u(maps.T_idx);
p  = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
pcn = repmat(p, maps.p_rep);

fdxcn = pr.ode_DFDX(pr, T0+T*mesh.tcn', xcn, pcn);
dxode = mesh.fdxka.*fdxcn;
dxode = sparse(maps.fdxrows, maps.fdxcols, dxode(:));
dxode = (0.5*T/maps.NTST)*dxode*maps.W-maps.Wp;

J0 = [dxode; maps.Q];

cdata = coco_get_chart_data(cseg.src_chart, tst.fid);
if ~isempty(cdata)
  M0 = [tst.row0; J0]\[cdata.M0; tst.rhs0];
else
  M0 = tst.M0;
end

M  = [tst.B1+M0'*tst.B2; J0]\tst.rhs;
It = 1;
while It<pr.coll.NBItMX && norm((M(:)-M0(:)))>0.1
  M0 = tst.la1*M0 + tst.la2*M;
  J  = [tst.B1+M0'*tst.B2; J0];
  M  = J\tst.rhs;
  It = It+1;
end
tst.M0  = tst.la1*M0 + tst.la2*M;
tst.row = tst.B1+M'*tst.B2;

data.coll_tst = tst;

end

function [prob, stat, xtr] = remesh(prob, data, chart, ub, Vb) %#ok<INUSD>

pr = data.pr;
seg  = pr.coll_seg;
tst  = pr.coll_tst;
dim  = seg.int.dim;
maps = seg.maps;
mesh = seg.mesh;

tst.row0    = [eye(dim) zeros(dim, maps.xbp_idx(end)-dim)];
tst.rhs0    = zeros(maps.xbp_idx(end)-dim, dim);
tst.M1_idx  = maps.xbp_idx(end-dim+(1:dim));
% tst.tbp     = mesh.tbp;
% tst.tbp_idx = maps.tbp_idx;
% tst.xbp_shp = maps.xbp_shp;

switch pr.coll.method
  case '2I'
    wts     = kron(mesh.gwt, eye(dim));
    kas     = kron(mesh.gka, eye(dim));
    tst.B1  = (0.5/maps.NTST)*((wts.*kas)*maps.W); % Check kappas!
    tst.rhs = [2*eye(dim); zeros(maps.xbp_idx(end)-dim, dim)];
    tst.la1 = 0.5;
    tst.la2 = 0.5;
  case '3I'
    ze = sparse(dim, maps.xbp_idx(end)-2*dim);
    tst.B1  = [speye(dim) ze speye(dim)];
    tst.rhs = [3*eye(dim); zeros(maps.xbp_idx(end)-dim, dim)];
    tst.la1 = 2/3;
    tst.la2 = 1-tst.la1; % don't use 1/3 [roundoff errors]
  otherwise
    error('%s: unknown mode', mfilename);
end
tst.B2 = (0.5/pr.coll.NTST)*maps.W'*mesh.wts2*mesh.kas2*maps.W; % Check kappas!

% u = coco_get_func_data(prob, seg.fid, 'u0');
% 
% x = u(maps.xbp_idx);
% T = u(maps.T_idx);
% p = u(maps.p_idx);
% 
% xx = reshape(maps.W*x, maps.x_shp);
% pp = repmat(p, maps.p_rep);
% 
% dxode = data.ode_DFDX(data, T*mesh.tcn', xx, pp);
% dxode = mesh.fdxka.*dxode;
% dxode = sparse(maps.fdxrows, maps.fdxcols, dxode(:));
% dxode = (0.5*T/maps.NTST)*dxode*maps.W-maps.Wp;
% 
% J0 = [dxode; maps.Q];
% 
% cdata = coco_get_chart_data(chart, tst.fid);
% M0 = [tst.row0; J0]\[cdata.M0; tst.rhs0];
% M  = [tst.B1+M0'*tst.B2; J0]\tst.rhs;
% It = 1;
% while It<data.coll.NBItMX && norm((M(:)-M0(:)))>0.1
%   M0 = tst.la1*M0 + tst.la2*M;
%   J  = [tst.B1+M0'*tst.B2; J0];
%   M  = J\tst.rhs;
%   It = It+1;
% end
% tst.M0  = tst.la1*M0 + tst.la2*M;
% tst.row = tst.B1+M'*tst.B2;

data.coll_tst = tst;

xtr  = [];
uidx = coco_get_func_data(prob, seg.fid, 'uidx');
prob = coco_change_func(prob, data, 'uidx', uidx);
stat = 'success';

end
