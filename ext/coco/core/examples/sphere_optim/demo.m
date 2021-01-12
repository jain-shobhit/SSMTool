%% Stationary points on a four-dimensional sphere
%
% This demo illustrates the successive continuation approach to searching
% for stationary points that satisfy the necessary conditions for a local
% extremum in the linear combination u1 + u2 + u3 + u4 along the
% three-dimensional sphere: u1^2 + u2^2 + u3^2 + u4^2 = 1.
%
% The augmented continuation problem constructed below has dimensional
% deficit -0-3=-3. The continuation parameters 'u2', 'u3', 'sum',
% 'd.u1', 'd.u2', 'd.u3', 'd.u4', and 'd.sum' are all initially inactive.
% Since each stage of continuation is along a one-dimensional manifold,
% four of the initially-inactive parameters must be released in each stage.
% In particular, 'sum' must always be active. Consequently, 'd.sum' should
% be active during the first two stages of continuation, the second of
% which terminates when this parameter equals 1. It should be inactive (and
% equal to 1) for all remaining stages of continuation. Since 'u1' and 'u4'
% are active during all stages of continuation, 'd.u1' and
% 'd.u4' should remain inactive during all stages of continuation. Since
% 'u2' is inactive during the first two stages of continuation and
% released only in the third stage of continuation, 'd.u2' should be active
% during all three stages and inactive thereafter. Similarly, since 'u3' is
% inactive during the first three stages of continuation and released only
% in the fourth stage of continuation, 'd.u3' should be active during all
% four stages.
%
% In the first stage of continuation, local extrema in 'sum' are located
% along a one-dimensional solution manifold with trivial Lagrange
% multipliers. Each of these is a branch point, from which emanates a
% secondary one-dimensional submanifold along which the Lagrange
% multipliers take on nontrivial values. As explained above, we terminate
% continuation along such a manifold when 'd.sum' equals 1. Stationary
% points within the computational domain correspond to points with
% vanishing 'd.u2' and 'd.u3'.

% The figure shows the one-dimensional solution manifolds obtained in the
% first, third, and fourth stages of continuation.

%% Initial encoding

prob = coco_prob;

fcns1 = { @sphere, @sphere_du, @sphere_dudu };
fcns2 = { @comb, @comb_du, @comb_dudu };

%% First run to find local extrema
% zero problem
prob1 = coco_add_func(prob, 'sphere', fcns1{:}, [], 'zero', ...
  'u0', [1;0;0;0]);

prob1 = coco_add_func(prob1, 'sum', fcns2{:}, [], 'inactive', ...
  'sum', 'uidx', (1:4)');

prob1 = coco_add_pars(prob1, 'pars1', [2; 3], {'u2', 'u3'}, 'inactive');

prob1 = coco_add_pars(prob1, 'pars2', [1; 4], {'u1', 'u4'}, 'active');

% adjoints
prob1 = coco_add_adjt(prob1, 'sphere');

prob1 = coco_add_adjt(prob1, 'sum', 'd.sum', 'aidx', (1:4)');

prob1 = coco_add_adjt(prob1, 'pars1', {'d.u2', 'd.u3'}, 'aidx', [2; 3]);

prob1 = coco_add_adjt(prob1, 'pars2', {'d.u1', 'd.u4'}, 'aidx', [1; 4]);
% continuation
cont_args = { 1, {'sum' 'd.sum' 'd.u2' 'd.u3'}, [1 2] };
bd1 = coco(prob1, 'sphere1', [], cont_args{:});

%% Switch at fold to branch with nontrivial multipliers
lab = coco_bd_labs(bd1, 'BP');

% branch switch data
chart = coco_read_solution('sphere1', lab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');
    
% zero problem
[chart, uidx] = coco_read_solution('sphere', 'sphere1', lab(1), ...
  'chart', 'uidx');

prob2 = coco_add_func(prob, 'sphere', fcns1{:}, [], 'zero', ...
  'u0', chart.x, 't0', cdata.v(uidx));

prob2 = coco_add_func(prob2, 'sum', fcns2{:}, [], 'inactive', ...
  'sum', 'uidx', (1:4)');

prob2 = coco_add_pars(prob2, 'pars1', [2; 3], {'u2', 'u3'}, 'inactive');

prob2 = coco_add_pars(prob2, 'pars2', [1; 4], {'u1', 'u4'}, 'active');

% adjoints
[chart, lidx] = coco_read_adjoint('sphere', 'sphere1', lab, ...
  'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'sphere', 'l0', chart.x, ...
  'tl0', cdata.v(lidx));
  
[chart, lidx] = coco_read_adjoint('sum', 'sphere1', lab, ...
  'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'sum', 'd.sum', 'aidx', (1:4)', ...
  'l0', chart.x, 'tl0', cdata.v(lidx));
  
dsum_int = [chart.x 1];
[chart, lidx] = coco_read_adjoint('pars1', 'sphere1', lab, ...
  'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'pars1', {'d.u2', 'd.u3'}, 'aidx', [2; 3], ...
  'l0', chart.x, 'tl0', cdata.v(lidx));
  
[chart, lidx] = coco_read_adjoint('pars2', 'sphere1', lab, ...
  'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'pars2', {'d.u1', 'd.u4'}, 'aidx', [1; 4], ...
  'l0', chart.x, 'tl0', cdata.v(lidx));
    
% continuation
cont_args = { 1, {'d.sum' 'sum' 'd.u2' 'd.u3'}, {[0 1], [-2 2]} };
bd2 = coco(prob2, 'sphere2', [], cont_args{:});

%% Constrain nontrivial multiplier and release additional continuation parameters
lab = coco_bd_labs(bd2, 'EP');
  
% zero problem
chart = coco_read_solution('sphere', 'sphere2', lab(2), 'chart');

prob3 = coco_add_func(prob, 'sphere', fcns1{:}, [], 'zero', ...
  'u0', chart.x);

prob3 = coco_add_func(prob3, 'sum', fcns2{:}, [], 'inactive', ...
  'sum', 'uidx', (1:4)');

prob3 = coco_add_pars(prob3, 'pars1', [2; 3], {'u2', 'u3'}, 'inactive');

prob3 = coco_add_pars(prob3, 'pars2', [1; 4], {'u1', 'u4'}, 'active');
    
% adjoints
chart = coco_read_adjoint('sphere', 'sphere2', lab(2), 'chart');
prob3 = coco_add_adjt(prob3, 'sphere', 'l0', chart.x);

chart = coco_read_adjoint('sum', 'sphere2', lab(2), 'chart');
prob3 = coco_add_adjt(prob3, 'sum', 'd.sum', 'aidx', (1:4)', ...
  'l0', chart.x);
  
chart = coco_read_adjoint('pars1', 'sphere2', lab(2), 'chart');
prob3 = coco_add_adjt(prob3, 'pars1', {'d.u2', 'd.u3'}, ...
  'aidx', [2; 3], 'l0', chart.x);
  
chart = coco_read_adjoint('pars2', 'sphere2', lab(2), 'chart');
prob3 = coco_add_adjt(prob3, 'pars2', {'d.u1', 'd.u4'}, ...
  'aidx', [1; 4], 'l0', chart.x);
  
% events
prob3 = coco_add_event(prob3, 'OPT', 'd.u2', 0);

% continuation
cont_args = { 1, {'d.u2', 'sum', 'u2', 'd.u3'}, {[] [-2 2]} };
bd3 = coco(prob3, 'sphere3', [], cont_args{:});

%% Constrain nontrivial adjoint and release additional continuation parameters
lab = coco_bd_labs(bd3, 'OPT');
    
% zero problem
chart = coco_read_solution('sphere', 'sphere3', lab, 'chart');

prob4 = coco_add_func(prob, 'sphere', fcns1{:}, [], 'zero', ...
  'u0', chart.x);
  
prob4 = coco_add_func(prob4, 'sum', fcns2{:}, [], 'inactive', ...
  'sum', 'uidx', (1:4)');
  
prob4 = coco_add_pars(prob4, 'pars1', [2; 3], {'u2', 'u3'}, 'inactive');

prob4 = coco_add_pars(prob4, 'pars2', [1; 4], {'u1', 'u4'}, 'active');
    
% adjoints
chart = coco_read_adjoint('sphere', 'sphere3', lab, 'chart');
prob4 = coco_add_adjt(prob4, 'sphere', 'l0', chart.x);
    
chart = coco_read_adjoint('sum', 'sphere3', lab, 'chart');
prob4 = coco_add_adjt(prob4, 'sum', 'd.sum', 'aidx', (1:4)', ...
  'l0', chart.x);

chart = coco_read_adjoint('pars1', 'sphere3', lab, 'chart');
prob4 = coco_add_adjt(prob4, 'pars1', {'d.u2', 'd.u3'}, ...
  'aidx', [2; 3], 'l0', chart.x);

chart = coco_read_adjoint('pars2', 'sphere3', lab, 'chart');
prob4 = coco_add_adjt(prob4, 'pars2', {'d.u1', 'd.u4'}, ...
  'aidx', [1; 4], 'l0', chart.x);

% events
prob4 = coco_add_event(prob4, 'OPT', 'd.u3', 0);

% continuation
cont_args = { 1, {'d.u3', 'sum', 'u2', 'u3'}, {[] [-2 2]} };
coco(prob4, 'sphere4', [], cont_args{:});

%% Graphical representation

theme = struct();
theme.lspec   = {{'g-', 'LineWidth', 1}, {'g-', 'LineWidth', 1}};
theme.special = {'OPT'};
theme.OPT     = {'kp', 'MarkerFaceColor', 'r', 'MarkerSize', 8};
      
figure(1)
clf
coco_plot_bd(theme, 'sphere1', 'u2', 'u3', 'u4')
hold on
coco_plot_bd(theme, 'sphere3', 'u2', 'u3', 'u4')
coco_plot_bd(theme, 'sphere4', 'u2', 'u3', 'u4')
hold off
grid on
view (-35, 40)
