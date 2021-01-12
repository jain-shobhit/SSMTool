function demo_arnold
% See http://dx.doi.org/10.1016/j.jcp.2006.05.041 for information about the
% golden-mean and Farey tree algorithms.

pq = tongues([1 4],[1 3],2);
at = pq(1,:)./pq(2,:);
[at idx] = sort(at);
pq = pq(:,idx);
h  = [0.2 1 0.3 0.7 0.3];
qp = sqp_curves(1/4,1/3,2);
[bd lab_at lab_qp] = run_po(at,qp);

for lab=lab_qp; run_qp(bd, lab); end
for lab=lab_at; run_at(bd, lab, at, pq, h); end

coco_use_recipes_toolbox

end

function [bd at qp] = run_po(at, qp)

if coco_exist({'AT' 'run_po'}, 'run')
  bd = coco_bd_read({'AT' 'run_po'});
else
  coco_use_recipes_toolbox coll_v1 po_v1
  
  p0 = [3.5;0.4;0];
  [t0 x0] = ode45(@(t,x) lang_red(x,p0), [0 5.3], [0.3;0.4]);
  
  TAT = 2*pi./(p0(1)*at);
  TQP = 2*pi./(p0(1)*qp);
  
  prob = coco_prob();
  prob = coco_set(prob, 'cont', 'PtMX', [15 0]);
  prob = coco_set(prob, 'coll', 'NTST', 20);
  prob = po_isol2orb(prob, '', @lang_red, t0, x0, {'om' 'ro' 'eps'}, p0);
  [data uidx] = coco_get_func_data(prob', 'po', 'data', 'uidx');
  prob = coco_add_pars(prob, '', uidx(data.x0_idx), {'x0' 'y0'}, 'active');
  prob = coco_add_event(prob, 'AT', 'po.period', TAT);
  prob = coco_add_event(prob, 'QP', 'po.period', TQP);
  bd = coco(prob, {'AT' 'run_po'}, [], 1, {'ro' 'po.period'}, {[0 1]});
end

lab = coco_bd_col(bd, 'LAB');
idx = coco_bd_idxs(bd, 'AT');
at  = fliplr([lab{idx}]);
idx = coco_bd_idxs(bd, 'QP');
qp  = fliplr([lab{idx}]);

end

function [rho]=sqp_curves(rho1, rho2, n)
% compute golden-mean subdivisions of [rho1, rho2] up to level n

rho = [rho1 rho2];
R   = rho;
gmr = 2/(1+sqrt(5));
gml = 1-gmr;

for i=1:n
  RR = [];
  m  = size(R,2);
  for j=1:m-1
    rho = [rho gmr*R(j)+gml*R(j+1) gml*R(j)+gmr*R(j+1)];
    RR = [RR R(j) gmr*R(j)+gml*R(j+1) gml*R(j)+gmr*R(j+1)];
  end
  R = [RR R(m)];
end

rho = rho(3:size(rho,2));

end

function [pq]=tongues(pq1, pq2, n)
% compute Farey-tree starting with pq1(1)/pq1(2) and pq2(1)/pq2(2) up to
% level n

P = [pq1(1) pq2(1)];
Q = [pq1(2) pq2(2)];

pq = [P; Q];

for i=1:n
  PP = []; QQ = [];
  m  = size(P,2);
  for j=1:m-1
    pq = [pq [P(j)+P(j+1); Q(j)+Q(j+1)]];
    PP = [PP P(j) P(j)+P(j+1)];
    QQ = [QQ Q(j) Q(j)+Q(j+1)];
  end
  P = [PP P(m)];
  Q = [QQ Q(m)];
end

end

function run_qp(bd, rlab)

runid = {'AT', sprintf('run_qp_%d',rlab)};
if ~coco_exist(runid, 'run')
  coco_use_recipes_toolbox coll_v1 msbvp_v1
  
  ro   = coco_bd_val(bd, rlab, 'ro');
  T_po = coco_bd_val(bd, rlab, 'po.period');
  x0   = coco_bd_val(bd, rlab, 'x0');
  y0   = coco_bd_val(bd, rlab, 'y0');
  
  p0     = [3.5;ro;0];
  N      = 75;
  tout   = linspace(0,T_po,2*N+2);
  [t x0] = ode45(@(t,x) lang_red(x,p0), tout, [x0;y0]);
  
  T_ret = 2*pi/p0(1);
  tt    = linspace(0,1,20*(2*N+1))';
  t1    = T_ret*tt;
  stt   = sin(tt*2*pi);
  ctt   = cos(tt*2*pi);
  coll_args = {};
  for i=1:2*N+1
    [t xx]  = ode45(@(t,x) lang_red(x,p0), [0 T_ret], x0(i,:));
    xx      = interp1(t, xx, t1);
    x1      = [ctt.*xx(:,1) stt.*xx(:,1) xx(:,2)];
    coll_args = [coll_args {@lang @lang_DFDX @langX_DFDP t1 x1 ...
      [p0; T_ret]}];
  end
  
  Th  = 2*pi*(0:2*N)/(2*N+1);
  Th  = kron(1:N, Th');
  F   = [ones(2*N+1,1) ...
    2*reshape([cos(Th); sin(Th)], [2*N+1 2*N])]'/(2*N+1);
  Th  = (1:N)*2*pi*T_ret/T_po;
  SIN = [zeros(size(Th)); sin(Th)];
  R   = diag([1, kron(cos(Th), [1, 1])])+diag(SIN(:),1)-diag(SIN(:),-1);
  data.F  = kron(F, eye(3));
  data.RF = kron(R*F, eye(3));
  
  prob = coco_prob();
  prob = coco_set(prob, 'coll', 'NTST', 20);
  prob = msbvp_isol2segs(prob, '', coll_args{:}, ...
    {'om' 'ro' 'eps' 'T_ret'}, @torus_bc, @torus_bc_DFDX, data);
  prob = coco_set(prob, 'corr', 'SubItMX', 2, 'ItMX', 30);
  coco(prob, {'AT','run0'}, [], 0, {'ro' 'T_ret'});
  
  prob = coco_prob();
  prob = coco_set(prob, 'coll', 'NTST', 20);
  prob = coco_set(prob, 'cont', 'h_max', 5);
  prob  = msbvp_sol2segs(prob, '', {'AT','run0'}, 1);
  prob = coco_set(prob, 'cont', 'ItMX', 70);
  coco(prob, runid, [], 1, {'eps' 'ro' 'T_ret'}, [-0.3 0.3]);
end

end


function y = lang_red(x, p)
%LANG_RED   'coll'-compatible encoding of symmetry-reduced langford vector field.

x1 = x(1,:);
x2 = x(2,:);
ro = p(2,:);

y(1,:) = (x2-0.7).*x1;
y(2,:) = 0.6+x2-x2.^3/3-x1.^2.*(1+ro.*x2);

end

function J = langX_DFDP(x, p)

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);

J = zeros(3,4,size(x,2));
J(1,1,:) = -x2;
J(2,1,:) = x1;
J(3,2,:) = -x3.*(x1.^2+x2.^2);
J(3,3,:) = x3.*x1.^3;

end

function y = torus_bc(data, T, x0, x1, p)
  y = [T-p(4); data.F*x1-data.RF*x0; x0(2); x0(4)-x0(1)];
end

function J = torus_bc_DFDX(data, T, x0, x1, p)

nt = numel(T);
nx = numel(x0);
np = numel(p);

J1 = zeros(2,nt+2*nx+np);
J1(1,nt+2) = 1;
J1(2,nt+[1 4]) = [-1 1];

J = [eye(nt), zeros(nt,2*nx+np-1), -ones(nt,1);
  zeros(nx,nt), -data.RF, data.F, zeros(nx,np);
  J1];

end

function run_at(bd, rlab, at, pq, h)
runid = {'AT', sprintf('run_at_%d',rlab)};
if ~coco_exist(runid, 'run')
  coco_use_recipes_toolbox atlas2d_v6 bvp_v1 coll_v1
  
  ro   = coco_bd_val(bd, rlab, 'ro');
  T_po = coco_bd_val(bd, rlab, 'po.period');
  x0   = coco_bd_val(bd, rlab, 'x0');
  y0   = coco_bd_val(bd, rlab, 'y0');
  
  p0 = [3.5;ro;0];
  rot = 2*pi/(p0(1)*T_po);
  [v idx] = min(abs(at-rot)); %#ok<ASGLU>
  
  [t x0] = ode45(@(t,x) lang(x,p0), [0 pq(1,idx)*T_po], [x0;0;y0]);
  
  prob = coco_prob();
  prob = coco_set(prob, 'coll', 'NTST', pq(2,idx)*20);
  prob = bvp_isol2seg(prob, '', @lang, @lang_DFDX, @lang_DFDP,  ...
    t, x0, {'om' 'ro' 'eps'}, p0, @po_bc, @po_bc_DFDX);
  
  prob = coco_add_slot(prob, 'AT_bddat', @add_AT, [], 'bddat');
  
  prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
  prob = coco_set(prob, 'cont', 'h', h(idx), 'almax', 30);
  prob = coco_set(prob, 'cont', 'PtMX', 1500);
  prob = coco_set(prob, 'cont', 'NPR', 300);
  coco(prob, runid, [], 2, {'ro' 'eps'}, {[] [-0.1 0.1]});
end

end

function [data res] = add_AT(prob, data, command, varargin)

res = {};
switch command
  case 'init'
    res   = 'AT';
  case 'data'
    chart = varargin{1};
    xidx = coco_get_func_data(prob, 'bvp.seg.coll', 'uidx');
    pidx = coco_get_func_data(prob, 'bvp.seg.coll.pars', 'uidx');
    res  = chart.x([pidx;xidx([1 2 3])]);
end

end
