%% Local minima in the presence of inequality constraints
%
% This demo illustrates the successive continuation approach to searching
% for local minima of a scalar-valued objective function on a region
% defined by a single scalar inequality.
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

% The figures show the one-dimensional solution manifolds obtained in the
% four stages of continuation.

%% Initial construction

addpath([pwd(),'/../../../contributed/symcoco/toolbox']);
u     = sym('u',[3,1]);
fu    = u(2)*(cos(u(3))+(1-u(1))*sin(u(3)))/((u(1)-1)^2+1);
sco_sym2funcs(fu, {u}, {'u'}, 'filename', 'sym_linode');
Fu    = sco_gen(@sym_linode);

fcn   = @(f) @(p,d,u) deal(d, f(u));
funcs = { fcn(Fu('')), fcn(Fu('u')), fcn(Fu({'u','u'})) };
bound = { fcn(@(u) u-1), fcn(@(u) 1), fcn(@(u) 0) };
data  = struct('funcs', {funcs}, 'bound', {bound});

prob  = coco_prob;
prob  = coco_set(prob, 'cont', 'PtMX', 50, 'NPR', inf);

%% First stage of continuation
probb = isol2prob(prob, data, [1; 2; 4]);
coco(probb, 'run1', [], 1, {'v' 'k' 'f' 'd.v' 'd.theta'});

%% Branch switch to secondary manifold
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

function prob = isol2prob(prob, data, u0)

prob = coco_add_func(prob, 'vel', data.funcs{:}, [], 'inactive', 'v', ...
  'u0', u0);
prob = coco_add_pars(prob, 'pars', 1:3, {'k' 'f' 'theta'});
prob = coco_add_func(prob, 'bound', data.bound{:}, [], ...
  'inequality', 'bd', 'uidx', 2);
prob = coco_add_adjt(prob, 'vel', 'd.v');
prob = coco_add_adjt(prob, 'pars', {'d.k' 'd.f' 'd.theta'}, 'aidx', 1:3);
prob = coco_add_adjt(prob, 'bound', 'ncp.f', 'aidx', 2);

end

function prob = BP2prob(prob, data, run, lab)

chart = coco_read_solution(run, lab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');

[chart, uidx] = coco_read_solution('vel', run, lab, 'chart', 'uidx');
prob = coco_add_func(prob, 'vel', data.funcs{:}, [], 'inactive', 'v', ...
  'u0', chart.x, 't0', cdata.v(uidx));
prob = coco_add_pars(prob, 'pars', 1:3, {'k' 'f' 'theta'});
prob = coco_add_func(prob, 'bound', data.bound{:}, [], ...
  'inequality', 'bd', 'uidx', 2);
[chart, lidx] = coco_read_adjoint('vel', run, lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'vel', 'd.v', 'l0', chart.x, ...
  'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('pars', run, lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'pars', {'d.k' 'd.f' 'd.theta'}, ...
  'aidx', 1:3, 'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('bound', run, lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'bound', 'ncp.f', 'aidx', 2, 'l0', chart.x, ...
  'tl0', cdata.v(lidx));

end

function prob = sol2prob(prob, data, run, lab)

chart = coco_read_solution('vel', run, lab, 'chart');
prob  = coco_add_func(prob, 'vel', data.funcs{:}, [], 'inactive', 'v', ...
  'u0', chart.x);
prob  = coco_add_pars(prob, 'pars', 1:3, {'k' 'f' 'theta'});
prob  = coco_add_func(prob, 'bound', data.bound{:}, [], ...
  'inequality', 'bd', 'uidx', 2);
chart = coco_read_adjoint('vel', run, lab, 'chart');
prob  = coco_add_adjt(prob, 'vel', 'd.v', 'l0', chart.x);
chart = coco_read_adjoint('pars', run, lab, 'chart');
prob  = coco_add_adjt(prob, 'pars', {'d.k' 'd.f' 'd.theta'}, ...
  'aidx', 1:3, 'l0', chart.x);
chart = coco_read_adjoint('bound', run, lab, 'chart');
prob  = coco_add_adjt(prob, 'bound', 'ncp.f', 'aidx', 2, 'l0', chart.x);

end
