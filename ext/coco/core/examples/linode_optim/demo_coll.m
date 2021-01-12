%% Local minima in the harmonically excited linear oscillator
%
% This demo illustrates the successive continuation approach to searching
% for local minima in the variable x2(0) along the three-dimensional
% manifold of periodic solutions to the system of differential equations
%
%     x1' = x2, x2' = -x2-k*x1+f*cos(t+theta)
%
% on the feasible region f<=1. 
%
% The augmented continuation problem constructed below has dimensional
% deficit -1-3=-4. By construction, 'v', 'k', 'f', 'theta', 'd.v', 'd.k',
% 'd.f', 'd.theta' are all initially inactive. Since each stage of
% continuation is along a one-dimensional manifold, five of these must be
% released in each stage. In particular, 'v' must always be active. Since
% 'k' and 'f' are released in all four stages, 'd.k' and 'd.f' remain fixed
% and equal to 0 throughout the analysis.
%
% In this version, the initial point lies in the infeasible region. In the
% first stage of continuation, local extrema in 'v' are located along a
% one-dimensional solution manifold with trivial Lagrange multipliers. Each
% of these is a branch point, from which emanates a secondary
% one-dimensional submanifold along which the Lagrange multipliers take on
% nontrivial values. As explained above, we terminate continuation along
% such a manifold when 'd.v' equals 1. Stationary points within the
% computational domain correspond to points with vanishing 'd.k', 'd.f',
% and 'd.theta'.

% The figure shows the one-dimensional solution manifolds obtained in the
% four stages of continuation.

%% Initial encoding

addpath([pwd(),'/../../../contributed/symcoco/toolbox']);
syms t x1 x2 k f theta
dxdt  = [x2; -x2-k*x1+f*cos(t+theta)];
sco_sym2funcs(dxdt, {t,[x1;x2],[k;f;theta]}, {'t','x','p'}, ...
  'vector', [0,1,1], 'filename', 'sym_linode');
F     = sco_gen(@sym_linode);

u     = sym('u',[6,1]); % x0(1:2), x1(1:2), T0, T
fu    = [u(3)-u(1); u(4)-u(2); u(5); u(6)-2*pi];
sco_sym2funcs(fu, {u}, {'u'}, 'filename', 'sym_bc');
Fu    = sco_gen(@sym_bc);

fcn   = @(f) @(p,d,u) deal(d, f(u));
funcs = { F(''), F('x'), F('p'), F('t'), ...
  F({'x','x'}), F({'x','p'}), F({'p','p'}), ...
  F({'t','x'}), F({'t','p'}), F({'t','t'})};
bcs   = { fcn(Fu('')), fcn(Fu('u')), fcn(Fu({'u','u'})) };
bound = { fcn(@(u) u-1), fcn(@(u) 1), fcn(@(u) 0) };
data  = struct('funcs', {funcs}, 'bcs', {bcs}, 'bound', {bound});

p0 = [1; 2; 4];
[t0, x0] = ode45(@(t,x) funcs{1}(t,x,p0), [0, 2*pi], [-1.5; -1.3]);

prob  = coco_prob;
prob  = coco_set(prob, 'ode', 'autonomous', false);
prob  = coco_set(prob, 'coll', 'NTST', 20);
prob  = coco_set(prob, 'cont', 'PtMX', 50, 'NPR', inf, 'NAdapt', 1);

%% First stage of continuation
probb = isol2prob(prob, data, t0, x0, p0);
coco(probb, 'run1', [], 1, {'v' 'k' 'f' 'd.v' 'd.theta'});

%% Switch at fold to branch with nontrivial multipliers
bd1   = coco_bd_read('run1');
BPlab = coco_bd_labs(bd1, 'BP');

probb = BP2prob(prob, data, 'run1', BPlab(1));
coco(probb, 'run2', [], 1, {'d.v' 'k' 'f' 'v' 'd.theta'}, [0 1]);

%% Fix d.v=1 and release theta, while driving d.theta to 0
probb = sol2prob(prob, data, 'run2', 2);
coco(probb, 'run3', [], 1, {'d.theta' 'k' 'f' 'theta' 'v'}, [-1 0]);

%% Fix d.theta=0 and release ncp.f and drive it to 0
probb = sol2prob(prob, data, 'run3', 4);
coco(probb, 'run4', [], 1, {'ncp.f' 'k' 'f' 'theta' 'v'}, [0 2]);

%% Graphical representation
figure(1); clf; hold on
thm = struct();
thm.special = {'BP'};
coco_plot_bd(thm, 'run1', 'k', 'theta', 'v')
coco_plot_bd(thm, 'run2', 'k', 'theta', 'v')
coco_plot_bd(thm, 'run3', 'k', 'theta', 'v')
coco_plot_bd(thm, 'run4', 'k', 'theta', 'v')
hold off; grid on; view(3)

figure(2); clf; hold on
thm = struct();
thm.special = {'BP'};
coco_plot_bd(thm, 'run1', 'k', 'f', 'v')
coco_plot_bd(thm, 'run2', 'k', 'f', 'v')
coco_plot_bd(thm, 'run3', 'k', 'f', 'v')
coco_plot_bd(thm, 'run4', 'k', 'f', 'v')
hold off; grid on; view(3)

rmpath([pwd(),'/../../../contributed/symcoco/toolbox']);

%% Composite constructors

function prob = isol2prob(prob, data, t0, x0, p0)

coll_args = [ data.funcs, {t0, x0, {'k' 'f' 'theta'}, p0} ];
prob = ode_isol2coll(prob, '', coll_args{:});
[fdata, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = fdata.coll_seg.maps;
prob = coco_add_func(prob, 'po', data.bcs{:}, [], 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx]));
prob = coco_add_func(prob, 'bound', data.bound{:}, [], ...
  'inequality', 'bd', 'uidx', uidx(maps.p_idx(2)));
prob = coco_add_pars(prob, 'vel', uidx(maps.x0_idx(2)), 'v');

prob = adjt_isol2coll(prob, '');
[fdata, axidx] = coco_get_adjt_data(prob, 'coll', 'data', 'axidx');
opt = fdata.coll_opt;
prob = coco_add_adjt(prob, 'po', 'aidx', ...
  axidx([opt.x0_idx; opt.x1_idx; opt.T0_idx; opt.T_idx]));
prob = coco_add_adjt(prob, 'bound', 'ncp.f', 'aidx', axidx(opt.p_idx(2)));
prob = coco_add_adjt(prob, 'vel', 'd.v', 'aidx', axidx(opt.x0_idx(2)));

end

function prob = BP2prob(prob, data, run, lab)

prob = ode_BP2coll(prob, '', run, lab);
[fdata, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = fdata.coll_seg.maps;
prob = coco_add_func(prob, 'po', data.bcs{:}, [], 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx]));
prob = coco_add_func(prob, 'bound', data.bound{:}, [], ...
  'inequality', 'bd', 'uidx', uidx(maps.p_idx(2)));
prob = coco_add_pars(prob, 'vel', uidx(maps.x0_idx(2)), 'v');

chart = coco_read_solution(run, lab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');

prob = adjt_BP2coll(prob, '', run, lab);
[chart, lidx] = coco_read_adjoint('po', run, lab, 'chart', 'lidx');
[fdata, axidx] = coco_get_adjt_data(prob, 'coll', 'data', 'axidx');
opt   = fdata.coll_opt;
prob   = coco_add_adjt(prob, 'po', 'aidx', ...
  axidx([opt.x0_idx; opt.x1_idx; opt.T0_idx; opt.T_idx]), ...
  'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('bound', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'bound', 'ncp.f', ...
  'aidx', axidx(opt.p_idx(2)), 'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('vel', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'vel', 'd.v', 'aidx', ...
  axidx(opt.x0_idx(2)), 'l0', chart.x, 'tl0', cdata.v(lidx));

end

function prob = sol2prob(prob, data, run, lab)

prob = ode_coll2coll(prob, '', run, lab);
[fdata, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
maps = fdata.coll_seg.maps;
prob = coco_add_func(prob, 'po', data.bcs{:}, [], 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.x1_idx; maps.T0_idx; maps.T_idx]));
prob = coco_add_func(prob, 'bound', data.bound{:}, [], ...
  'inequality', 'bd', 'uidx', uidx(maps.p_idx(2)));
prob = coco_add_pars(prob, 'vel', uidx(maps.x0_idx(2)), 'v');

prob = adjt_coll2coll(prob, '', run, lab);
chart = coco_read_adjoint('po', run, lab, 'chart');
[fdata, axidx] = coco_get_adjt_data(prob, 'coll', 'data', 'axidx');
opt = fdata.coll_opt;
prob   = coco_add_adjt(prob, 'po', 'aidx', ...
  axidx([opt.x0_idx; opt.x1_idx; opt.T0_idx; opt.T_idx]), 'l0', chart.x);
chart = coco_read_adjoint('bound', run, lab, 'chart');
prob = coco_add_adjt(prob, 'bound', 'ncp.f', ...
  'aidx', axidx(opt.p_idx(2)), 'l0', chart.x);
chart = coco_read_adjoint('vel', run, lab, 'chart');
prob = coco_add_adjt(prob, 'vel', 'd.v', 'aidx', ...
  axidx(opt.x0_idx(2)), 'l0', chart.x);

end
