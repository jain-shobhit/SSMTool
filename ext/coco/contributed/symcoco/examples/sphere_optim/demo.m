%% Test symcoco for demo of stationary points on a four-dimensional sphere
%
% (Copied from |coco_folder/core/examples/sphere_optim|.) See that folder
% and CORE-tutorial for details about the demo. This file shows in its
% first parts how one can use symbolic derivatives of up to second order
% for coco computations if the expected function format is
% 
%   function [data,y_or_J]=func(prob,data,u)
% 
% After the first three parts the demo is identical to the original CORE
% demo |sphere_optim|.
%% Load path to symcoco routines
% (only |sco_gen| needed here)
clear
addpath([pwd(),'/../../toolbox']);
%% Wrapper converting simple functions to coco format
% The construction for |fcn| below takes a function |f| with a single input
% |u| and output |f(u)| and converts it into a function that has inputs
% |p|, |d| and |u| and returns |d| and |f(u)|. |obj| and |con| are
% shortcuts for creating right-hand sides and derivatives of the constraint
% function, defining the sphere (|con|) and the objective functional
% (|obj|). See <gen_sym_sphere.html>.
fcn = @(f) @(p,d,u) deal(d, f(u));
obj=sco_gen(@sym_sphere_obj);
con=sco_gen(@sym_sphere_constraint);

%% Initial encoding
% Note how the function and partial derivatives of |con| are used directly
% in line 3 below, as they are treated as r.h.s. of an equilibrium problem.
% In contrast, in line 4 below the objective functional |obj| is in core
% |coco| format, such that one has to apply the wrapper |fcn|, defined
% above. Apart from these two lines, the remainder of the demo is identical
% to the core |coco| demo |sphere_optim|.
prob = coco_prob;
prob = coco_set(prob, 'ode', 'vectorized', false);
funcs1 = { con(''),con('x'),con('p')};
funcs2 = { fcn(obj('')),fcn(obj('u')),fcn(obj({'u','u'}))};

%% First run to find local extrema
% zero problem
prob1 = ode_isol2ep(prob, '', funcs1{:}, 1, {'u2' 'u3' 'u4'}, [0; 0; 0]);

prob1 = coco_add_func(prob1, 'sum', funcs2{:}, [], ...
  'inactive', 'sum', 'uidx', 1:4);

prob1 = coco_add_pars(prob1, 'pars', 1, 'u1', 'active');

% adjoints
prob1 = adjt_isol2ep(prob1, '');

prob1 = coco_add_adjt(prob1, 'sum', 'd.sum', 'aidx', 1:4);

prob1 = coco_add_adjt(prob1, 'pars', 'd.u1', 'aidx', 1);

% continuation
cont_args = {1, {'sum' 'd.sum' 'd.u2' 'd.u3' 'u4'}, [1 2]};
bd1 = coco(prob1, 'sphere1', [], cont_args{:});

%% Switch at fold to branch with nontrivial multipliers
BPlab = coco_bd_labs(bd1, 'BP');

% zero problem
prob2 = ode_BP2ep(prob, '', 'sphere1', BPlab(1));

prob2 = coco_add_func(prob2, 'sum', funcs2{:}, [], ...
  'inactive', 'sum', 'uidx', 1:4);

prob2 = coco_add_pars(prob2, 'pars', 1, 'u1', 'active');

% branch switch data
chart = coco_read_solution('sphere1', BPlab(1), 'chart');
cdata = coco_get_chart_data(chart, 'lsol');

% adjoints
prob2 = adjt_BP2ep(prob2, '', 'sphere1', BPlab(1));
[chart, lidx] = coco_read_adjoint('sum', 'sphere1', BPlab(1), ...
  'chart', 'lidx');

prob2 = coco_add_adjt(prob2, 'sum', 'd.sum', 'aidx', 1:4, ...
  'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('pars', 'sphere1', BPlab(1), ...
  'chart', 'lidx');

prob2 = coco_add_adjt(prob2, 'pars', 'd.u1', 'aidx', 1, ...
  'l0', chart.x, 'tl0', cdata.v(lidx));

% continuation
cont_args = {1, {'d.sum' 'sum' 'd.u2' 'd.u3' 'u4'}, {[0 1], [-2 2]}};
bd2 = coco(prob2, 'sphere2', [], cont_args{:});

%% Constrain nontrivial multiplier and release additional continuation parameters
lab = coco_bd_labs(bd2, 'EP');

% zero problem
prob3 = ode_ep2ep(prob, '', 'sphere2', lab(2));

prob3 = coco_add_func(prob3, 'sum', funcs2{:}, [], ...
  'inactive', 'sum', 'uidx', 1:4);

prob3 = coco_add_pars(prob3, 'pars', 1, 'u1', 'active');

% adjoints
prob3 = adjt_ep2ep(prob3, '', 'sphere2', lab(2));

chart = coco_read_adjoint('sum', 'sphere2', lab(2), 'chart');
prob3 = coco_add_adjt(prob3, 'sum', 'd.sum', 'aidx', 1:4, 'l0', chart.x);

chart = coco_read_adjoint('pars', 'sphere2', lab(2), 'chart');
prob3 = coco_add_adjt(prob3, 'pars', 'd.u1', 'aidx', 1, 'l0', chart.x);

% events
prob3 = coco_add_event(prob3, 'OPT', 'd.u2', 0);

% continuation
cont_args = {1, {'d.u2' 'sum' 'u2' 'd.u3' 'u4'}, {[], [-2 2]}};
bd3 = coco(prob3, 'sphere3', [], cont_args{:});

%% constrain nontrivial adjoint and release additional continuation parameters
lab = coco_bd_labs(bd3, 'OPT');

% zero problem
prob4 = ode_ep2ep(prob, '', 'sphere3', lab(1));

prob4 = coco_add_func(prob4, 'sum', funcs2{:}, [], ...
  'inactive', 'sum', 'uidx', 1:4);

prob4 = coco_add_pars(prob4, 'pars', 1, 'u1', 'active');

% adjoints
prob4 = adjt_ep2ep(prob4, '', 'sphere3', lab(1));

chart = coco_read_adjoint('sum', 'sphere3', lab(1), 'chart');
prob4 = coco_add_adjt(prob4, 'sum', 'd.sum', 'aidx', 1:4, 'l0', chart.x);

chart = coco_read_adjoint('pars', 'sphere3', lab(1), 'chart');
prob4 = coco_add_adjt(prob4, 'pars', 'd.u1', 'aidx', 1, 'l0', chart.x);

% events
prob4 = coco_add_event(prob4, 'OPT', 'd.u3', 0);

% continuation
cont_args = {1, {'d.u3' 'sum' 'u2' 'u3' 'u4'}, {[], [-2 2]}};
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
%% Remove path to |symcoco|
rmpath([pwd(),'/../../toolbox']);
