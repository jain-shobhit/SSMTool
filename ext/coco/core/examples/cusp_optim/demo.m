%% Stationary points in the cusp normal form
%
% This demo illustrates the successive continuation approach to searching
% for stationary points that satisfy the necessary conditions for a local
% extremum in the variable u2 along the two-dimensional surface 
% u2 = u1*(u3-u1^2).
%
% The augmented continuation problem constructed below has dimensional
% deficit 0-2=-2. The continuation parameters 'u2', 'u3', 'd.u2', and
% 'd.u3' are all initially inactive. Since each stage of continuation is
% along a one-dimensional manifold, three of these parameters must be
% released in each stage. In particular, 'u2' must always be active.
% Consequently, 'd.u2' should be active during the first two stages of
% continuation, the second of which terminates when this parameter equals
% 1. It should be inactive (and equal to 1) for all remaining stages of
% continuation. Since 'u3' is inactive during the first two stages of
% continuation and released only in the third stage of continuation, 'd.u3'
% should be active during all three stages.
%
% In the first stage of continuation, local extrema in 'u2' are located
% along a one-dimensional solution manifold with trivial Lagrange
% multipliers. Each of these is a branch point, from which emanates a
% secondary one-dimensional submanifold along which the Lagrange
% multipliers take on nontrivial values. As explained above, we terminate
% continuation along such a manifold when 'd.u2' equals 1. Stationary
% points within the computational domain correspond to points with
% vanishing 'd.u3'.
%
% In this example, the third stage of continuation is along a family of
% fold points with a cusp coincident with the stationary point.

% The figure shows the one-dimensional solution manifolds obtained in the
% three stages of continuation. The folds along the first manifold are not
% visible in this projection.

%% Initial encoding
prob = coco_prob;

%% First run to find local extrema
% zero problem
prob1 = coco_add_func(prob, 'cusp', @cusp, @cusp_du, @cusp_dudu, [], ...
  'zero', 'u0', [0;0;0.5]);
prob1 = coco_add_pars(prob1, 'pars', [2; 3], {'u2', 'u3'}, 'inactive');

% adjoints
prob1 = coco_add_adjt(prob1, 'cusp');
prob1 = coco_add_adjt(prob1, 'pars', {'d.u2', 'd.u3'}, 'aidx', [2; 3]);

% continuation
bd1 = coco(prob1, 'cusp1', [], 1, {'u2', 'd.u2', 'd.u3'}, [-0.5 0.5]);

%% Switch at fold to branch with nontrivial multipliers
BPlab = coco_bd_labs(bd1, 'BP');

% zero problem
[chart, uidx] = coco_read_solution('cusp', 'cusp1', BPlab(1), ...
  'chart', 'uidx');
cdata = coco_get_chart_data(chart, 'lsol');

prob2 = coco_add_func(prob, 'cusp', @cusp, @cusp_du, @cusp_dudu, [], ...
  'zero', 'u0', chart.x, 't0', cdata.v(uidx));
prob2 = coco_add_pars(prob2, 'pars', [2; 3], {'u2', 'u3'}, 'inactive');

% adjoints
[chart, lidx] = coco_read_adjoint('cusp', 'cusp1', BPlab(1), ...
  'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'cusp', 'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('pars', 'cusp1', BPlab(1), ...
  'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'pars', {'d.u2', 'd.u3'}, 'aidx', [2; 3], ...
  'l0', chart.x, 'tl0', cdata.v(lidx));

% continuation
cont_args = {1, {'d.u2', 'u2', 'd.u3'}, {[0 1], [-0.5 0.5]}};
bd2 = coco(prob2, 'cusp2', [], cont_args{:});

%% Constrain nontrivial multiplier and release additional continuation parameters
lab = coco_bd_labs(bd2, 'EP');

% zero problem
chart = coco_read_solution('cusp', 'cusp2', lab(2), 'chart');
prob3 = coco_add_func(prob, 'cusp', @cusp, @cusp_du, @cusp_dudu, [], ...
  'zero', 'u0', chart.x);
prob3 = coco_add_pars(prob3, 'pars', [2; 3], {'u2', 'u3'}, 'inactive');

% adjoints
chart = coco_read_adjoint('cusp', 'cusp2', lab(2), 'chart');
prob3 = coco_add_adjt(prob3, 'cusp', 'l0', chart.x);
chart = coco_read_adjoint('pars', 'cusp2', lab(2), 'chart');
prob3 = coco_add_adjt(prob3, 'pars', {'d.u2', 'd.u3'}, 'aidx', [2;3], ...
  'l0', chart.x);

% events
prob3 = coco_add_event(prob3, 'OPT', 'd.u3', 0);

% continuation
cont_args = {1, {'d.u3', 'u2', 'u3'}, {[], [-0.5 0.5], [-2 2]}};
coco(prob3, 'cusp3', [], cont_args{:});

%% Graphical representation

figure(1); clf; hold on
coco_plot_bd('cusp1', 'u2', 'd.u3', 'd.u2')
coco_plot_bd('cusp2', 'u2', 'd.u3', 'd.u2')
thm = struct();
thm.lspec = {'g-', 'LineWidth', 1};
thm.special = {'OPT'};
thm.OPT = {'kp', 'MarkerFaceColor', 'r', 'MarkerSize', 8};
coco_plot_bd(thm, 'cusp3', 'u2', 'd.u3', 'd.u2')
hold off; grid on; view(-70,45)
