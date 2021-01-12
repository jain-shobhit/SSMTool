%% Catenary problem from Section 7.3.1 of Recipes for Continuation
%
% A technique for constructing initial solutions to nonlinear
% boundary-value problems is continuation in the interval length T. Here,
% one starts with a short trajectory segment for small T<<1 and continues
% in T until reaching the desired value. This process is called growing an
% initial orbit and is illustrated in Sect. 7.3.1 of Recipes for
% Continuation in the context of the problem from calculus of variations,
% introduced in Chap. 1.

% Figure 1 shows an initial solution guess obtained with an Euler step at
% t=0 with step size h=0.04, or by single point, and compares this
% approximation with the exact solution. The grown solution is shown in
% Figure 2. From this solution one can, finally, start a continuation in
% y1e. Figure 3 shows the labeled solutions obtained in this final
% continuation run with the trajectory with minimum value of 'y1e'
% highlighted in red.

%% Encoding

% The continuation problem encoded below includes four monitor functions
% that evaluate to the components of the trajectory end point at t=0, the
% first component of the trajectory end point at t=1, and the interval
% length, respectively, and four corresponding inactive continuation
% parameters 'y1s', 'y2s', 'y1e', and 'T'. Its dimensional deficit equals
% -1. The call to the coco entry-point function indicates a desired
% manifold dimension of 1. To this end, the continuation parameter 'T' and
% 'y1' are released and allowed to vary during continuation.

% Vectorized vector field, optional Jacobians are omitted.
cat = @(x,p) [x(2,:); (1+x(2,:).^2)./x(1,:)];

% Initial solution guess is given in Recipes for Continuation by single
% Euler step of duration 0.04.

% t0 = [0 0.04];
% x0 = [1 0; 1 0.04];

% As an alternative, initial solution guess can consist of a single point.
t0 = 0;
x0 = [1 0];

prob = coco_prob();
prob = ode_isol2coll(prob, '', cat, t0, x0, []);
data = coco_get_func_data(prob, 'coll', 'data'); % Extract function data
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'pars', ...
  [maps.x0_idx; maps.x1_idx(1); maps.T_idx], ...
  {'y1s' 'y2s' 'y1e' 'T'});
prob = coco_set(prob, 'cont', 'NAdapt', 10); % Remesh every ten steps
cont_args = {1, {'T' 'y1e' 'coll.err' 'coll.err_TF'}, [t0(end) 1]};

fprintf('\n Run=''%s'': Continue trajectory segments until T=1.\n', ...
  'coll1');

bd = coco(prob, 'coll1', [], cont_args{:});

%% Restarting continuation along different submanifold

% The continuation problem encoded below is identical to that above, but
% constructed from a stored solution from the previous run. The call to the
% coco entry-point function indicates a desired manifold dimension of 1. To
% this end, the continuation parameters 'y1e' and 'y2s' are released and
% allowed to vary during continuation.

labs = coco_bd_labs(bd, 'EP');

prob = coco_prob();
prob = ode_coll2coll(prob, '', 'coll1', labs(end));
data = coco_get_func_data(prob, 'coll', 'data');
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'pars', ...
  [maps.x0_idx; maps.x1_idx(1); maps.T_idx], ...
  {'y1s' 'y2s' 'y1e' 'T'});
prob = coco_set(prob, 'cont', 'NAdapt', 10, 'PtMX', 200);
cont_args = {1, {'y1e' 'y2s' 'coll.err' 'coll.err_TF'}, [0 3]};

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s''.\n', ...
  'coll2', labs(end), 'coll1');

coco(prob, 'coll2', [], cont_args{:});

%% Graphical representation of stored solutions

% Exact solution and initial guess
figure(1); clf; hold on; grid on; box on; axis([0 1 0 2])
t = linspace(0,1,1000);
x = cosh(t);
plot(t, x, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
plot(t0, x0(:,1), 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)
xlabel('t')
ylabel('x_1')
hold off

% Initial result of shooting method
figure(2); clf; hold on; grid on; box on; axis([0 1 0 2])
coco_plot_sol('coll1', labs(end), '')
hold off

% "One-parameter continuation" along family of trajectories
figure(3); clf; hold on; grid on; box on; axis([0 1 0 3])
thm = struct('special', {{'FP'}});
thm.sol.FP = {'r-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12};
coco_plot_sol(thm, 'coll2',  '')
hold off
