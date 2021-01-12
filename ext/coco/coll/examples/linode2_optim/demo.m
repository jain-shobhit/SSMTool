%% Stationary points in the harmonically excited linear oscillator
%
% This demo illustrates the successive continuation approach to searching
% for stationary points that satisfy the necessary conditions for a local
% extremum in the variable x2(0) along the two-dimensional manifold of
% solutions to the system of boundary-value problem
%
%    x1' = x2, x2' = -x2-k*x1+cos(x3+theta), x3' = 1
%
% where x1(2*pi)-x1(0)=0, x2(2*pi)-x2(0)=0, x3(2*pi)-x3(0)-2*pi=0, x3(0)=0.
%
% The augmented continuation problem constructed below has dimensional
% deficit -1-2=-3. The continuation parameters 'k', 'th', 'v', 'd.k',
% 'd.th', and 'd.v' are all initially inactive. Since each stage of
% continuation is along a one-dimensional manifold, four of these
% parameters must be released in each stage. In particular, 'v' must always
% be active. Consequently, 'd.v' should be active during the first two
% stages of continuation, the second of which terminates when this
% parameter equals 1. It should be inactive (and equal to 1) for all
% remaining stages of continuation. Since 'k' is active during all stages
% of continuation, 'd.k' should be inactive during all three stages. Since
% 'th' is inactive during the first two stages of continuation and released
% only in the third stage of continuation, 'd.th' should be active during
% all three stages.
%
% In the first stage of continuation, a local extremum in 'v' is located
% along a one-dimensional solution manifold with trivial Lagrange
% multipliers. This is a branch point, from which emanates a secondary
% one-dimensional submanifold along which the Lagrange multipliers take on
% nontrivial values. As explained above, we terminate continuation along
% this manifold when 'd.v' equals 1. Stationary points within the
% computational domain correspond to points with vanishing 'd.k' and
% 'd.th'.

% The figure shows the one-dimensional solution manifolds obtained in the
% first and third stages of continuation.

%% Initial encoding

prob = coco_prob;

%% First run to find local extrema
prob1 = coco_set(prob, 'coll', 'NTST', 25);

% zero problems
[t0, x0]  = ode45(@(t,x) linode(x, [0.98; 0.3]), [0 2*pi], ...
  [0.276303; 0.960863; 0]);
coll_args = {@linode, @linode_dx, @linode_dp, @linode_dxdx, ...
  @linode_dxdp, @linode_dpdp, t0, x0, {'k' 'th'}, [0.98; 0.3]};
prob1 = ode_isol2coll(prob1, '', coll_args{:});

[data, uidx] = coco_get_func_data(prob1, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
bc_funcs = {@linode_bc, @linode_bc_du, @linode_bc_dudu};
prob1 = coco_add_func(prob1, 'po', bc_funcs{:}, [], 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.x1_idx]));

prob1 = coco_add_pars(prob1, 'vel', uidx(maps.x0_idx(2)), 'v');

% adjoints
prob1 = adjt_isol2coll(prob1, '');

[data, axidx] = coco_get_adjt_data(prob1, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob1 = coco_add_adjt(prob1, 'po', 'aidx', ...
  axidx([opt.x0_idx; opt.x1_idx]));

prob1 = coco_add_adjt(prob1, 'vel', 'd.v', 'aidx', axidx(opt.x0_idx(2)));

% continuation
cont_args = {1, {'v' 'k' 'd.v' 'd.th'}, [0.9 2]};
coco(prob1, 'linode1', [], cont_args{:});

%% Switch at fold to branch with nontrivial multipliers
bd1   = coco_bd_read('linode1');
BPlab = coco_bd_labs(bd1, 'BP');

% zero problems
prob2 = ode_BP2coll(prob, '', 'linode1', BPlab(1));

[data, uidx] = coco_get_func_data(prob2, 'coll', 'data', 'uidx');
maps  = data.coll_seg.maps;
prob2 = coco_add_func(prob2, 'po', bc_funcs{:}, data, 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.x1_idx]));

prob2 = coco_add_pars(prob2, 'vel', uidx(maps.x0_idx(2)), 'v');

% branch switch data
chart = coco_read_solution('linode1', BPlab(1), 'chart');
cdata = coco_get_chart_data(chart, 'lsol');

% adjoints
prob2 = adjt_BP2coll(prob2, '', 'linode1', BPlab(1));

[chart, lidx] = coco_read_adjoint('po', 'linode1', BPlab(1), ...
  'chart', 'lidx');
[data, axidx] = coco_get_adjt_data(prob2, 'coll', 'data', 'axidx');
opt   = data.coll_opt;
prob2   = coco_add_adjt(prob2, 'po', 'aidx', ...
  axidx([opt.x0_idx; opt.x1_idx]), 'l0', chart.x, 'tl0', cdata.v(lidx));

[chart, lidx] = coco_read_adjoint('vel', 'linode1', BPlab(1), ...
  'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'vel', 'd.v', 'aidx', ...
  axidx(opt.x0_idx(2)), 'l0', chart.x, 'tl0', cdata.v(lidx));

% continuation
cont_args = {1, {'d.v' 'v' 'k' 'd.th'}, {[0 1], [.9 2]}};
coco(prob2, 'linode2', [], cont_args{:});

%% Constrain nontrivial multiplier and release additional continuation parameters
bd2 = coco_bd_read('linode2');
lab = coco_bd_labs(bd2, 'EP');

% zero problem
prob3 = ode_coll2coll(prob, '', 'linode2', lab(2));

[data, uidx] = coco_get_func_data(prob3, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob3 = coco_add_func(prob3, 'po', bc_funcs{:}, data, 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.x1_idx]));

prob3 = coco_add_pars(prob3, 'vel', uidx(maps.x0_idx(2)), 'v');

% adjoint
prob3 = adjt_coll2coll(prob3, '', 'linode2', lab(2));

chart = coco_read_adjoint('po', 'linode2', lab(2), 'chart');
[data, axidx] = coco_get_adjt_data(prob3, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob3   = coco_add_adjt(prob3, 'po', 'aidx', ...
  axidx([opt.x0_idx; opt.x1_idx]), 'l0', chart.x);

chart = coco_read_adjoint('vel', 'linode2', lab(2), 'chart');
prob3 = coco_add_adjt(prob3, 'vel', 'd.v', 'aidx', ...
  axidx(opt.x0_idx(2)), 'l0', chart.x);

% events
prob3 = coco_add_event(prob3, 'OPT', 'd.th', 0);

% continuation
cont_args = {1, {'d.th' 'v' 'k' 'th'}, {[], [.9 2]}};
coco(prob3, 'linode3', [], cont_args{:});

%% Graphical representation

figure(1); clf; hold on
thm = struct();
thm.special = {'BP'};
coco_plot_bd(thm, 'linode1', 'k', 'th', 'v')
thm.special = {'OPT'};
thm.OPT     = {'kp', 'MarkerFaceColor', 'r', 'MarkerSize', 8};
coco_plot_bd(thm, 'linode3', 'k', 'th', 'v')
hold off; grid on; view(3)
