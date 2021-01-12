%% Piecewise-smooth dynamical system from Example 9.2 in Recipes for Continuation
%
% The vector fields
%
% F(x,p;left)  = (-x2+(1-r)*x1, x1+(1-r)*x2)
% F(x,p;right) = (al*(be-r)*x1-(ga+be-r)*x2, al*(be-r)*x2+(ga+be-r)*x1)
%
% are associated with trajectory segments that may be stitched together to
% form closed curves in a planar state space. Here, we consider stitching
% of trajectory segments across the x1=0 axis, corresponding to periodic
% orbits in the piecewise-smooth dynamical system with F(x,p;left) applied
% when x1<0 and F(x,p;right) applied when x1>0.

% Panel (a) shows a family of two-segment periodic orbits obtained under
% variation in be. The dashed curve represents the initial solution guess
% and the gray curve represents the first point located on the solution
% manifold. Panel (b) shows a family of two-segment periodic orbits
% obtained under variation in al by restarting continuation from one of the
% orbits located in the first run.

%% Initial encoding

% The continuation problem structure encoded below includes three monitor
% functions that evaluate to the problem parameters, and the corresponding
% inactive continuation parameters 'al', 'be', and 'ga'. Its dimensional
% deficit equals 0. To obtain a one-dimensional family of periodic orbits,
% the continuation parameter 'be' is released and allowed to vary during
% continuation.

% Construct the initial solution guess. The segment signature, contained in
% the modes, events, and resets variables uniquely defines the multisegment
% boundary-value problem. Continuation is restricted to two-segment
% periodic orbits with switching surface along the vertical axis.

p0     = [1; 2; 1.5];
modes  = {'left' 'right'};
events = {'boundary' 'boundary'};
resets = {'switch' 'switch'};
t0 = linspace(0, pi, 100)';
x1 = [-sin(t0) 0.5+cos(t0)];
x2 = [ sin(t0) 0.5-cos(t0)];
t0 = {t0 0.5*t0};
x0 = {x1 x2};

% Anonymous, in-line definition of event and reset functions and their
% Jacobians. There is no value in vectorizing these, since they are never
% called with multiple points.

stop    = @(x,p,e) x(1);
stop_dx = @(x,p,e) [1 0];
stop_dp = @(x,p,e) zeros(1,3);
jump    = @(x,p,r) x;
jump_dx = @(x,p,r) eye(2);
jump_dp = @(x,p,r) zeros(2,3);

% Use segment-specific settings for the initial values of the
% discretization parameters NTST and NCOL of the COLL toolbox.

prob = coco_prob();
prob = coco_set(prob, 'hspo.orb.bvp.seg1.coll', 'NTST', 10, 'NCOL', 6);
prob = coco_set(prob, 'hspo.orb.bvp.seg2.coll', 'NTST', 20, 'NCOL', 4);
prob = ode_isol2hspo(prob, '', ...
  {@piecewise, stop, jump}, ...
  {@piecewise_dx, stop_dx, jump_dx}, ...
  {@piecewise_dp, stop_dp, jump_dp}, ...
  modes, events, resets, t0, x0, {'al' 'be' 'ga'}, p0);

% Allow adaptive remeshing of the trajectory discretization every 5 steps.
prob = coco_set(prob, 'cont', 'NAdapt', 5);

fprintf('\n Run=''%s'': Continue family of two-segment periodic orbits.\n', ...
  'pw1');

coco(prob, 'pw1', [], 1, 'be', [0 5]);


%% Restarting continuation

% The continuation problem encoded below is identical to that above, but
% constructed from stored data obtained in the previous run. Toolbox
% settings encoded in the stored toolbox data override the default values.
% A one-dimensional family of periodic orbits is obtained by releasing 'al'
% and allowing it to vary during continuation.

prob = ode_hspo2hspo(coco_prob(), '', 'pw1', 9);
prob = coco_set(prob, 'cont', 'NAdapt', 5);

fprintf(...
  '\n Run=''%s'': Continue two-segment periodic orbits from point %d in run ''%s''.\n', ...
  'pw2', 9, 'pw1');

coco(prob, 'pw2', [], 1, 'al', [0 4]);

%% Graphical representation of result

% Plot data: panel (a)
figure(1); clf; hold on

plot([0 0], [-1.5 5], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', [0.5 0.5 0.5])
coco_plot_sol('pw1', '', 'x', 'x')
thm = struct();
thm.sol.RO = {'Color', [0.5 0.5 0.5], 'LineWidth', 1};
coco_plot_sol(thm, 'pw1', 9, '', 'x', 'x')
plot(x1(:,1), x1(:,2), 'k', 'LineStyle', '--', 'LineWidth', 2)
plot(x2(:,1), x2(:,2), 'k', 'LineStyle', '--', 'LineWidth', 2)

hold off; grid on; box on; axis([-1.5 3.5 -1.5 5])

% Plot data: panel (b)
figure(2); clf; hold on

plot([0 0], [-4 5], 'LineStyle', '-', 'LineWidth', 1, ...
  'Color', [0.3 0.3 0.3])
coco_plot_sol('pw2', '', 'x', 'x')
thm = struct();
thm.sol.RO = {'Color', [0.5 0.5 0.5], 'LineWidth', 1};
coco_plot_sol(thm, 'pw2', 6, '', 'x', 'x')

hold off; grid on; box on; axis([-1.5 5 -4 5])
