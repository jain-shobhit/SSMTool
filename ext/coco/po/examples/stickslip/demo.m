%% An stick-slip oscillator from Examples 9.3 and 9.4 in Recipes for Continuation
%
% We investigate families of periodic orbits in a two-degree-of-freedom
% mechanical system in the presence of nonlinear parametric excitation,
% velocity jumps associated with instantaneous impacts with a mechanical
% constraint and transitions between sticking and slipping motion of one of
% the degrees of freedom in the presence of dry friction. The example
% illustrates a multi-segment boundary-value problem with trajectory
% segments governed by vector fields with different state-space dimensions.

% Panels (a), (b), and (c) shows solutions to a two-segment periodic orbit
% continuation problem under variations in the excitation amplitude V, but
% the solution in (c) is not physically realizable, since it includes
% excursion beyond the mechanical constraint. Panels (d), (e), and (f) show
% solutions to a three-segment periodic orbit continuation problem under
% variations in V. The solution in (f) is not physically realizable, since
% it includes excursion beyond the mechanical constraint. A short sliding
% segment is seen in panel (d). This shrinks to a single point, with zero
% duration, in the grazing periodic orbit shown in panel (e).

%% Initial encoding

% The continuation problem structure encoded below includes five monitor
% functions that evaluate to the problem parameters, and the corresponding
% inactive continuation parameters 'V', 'c', 'n', 'w', and 'e'. Its
% dimensional deficit equals 0. A one-dimensional family of non-impacting
% sticking periodic orbits is obtained by releasing 'V' and allowing this
% to vary during continuation.

% Construct initial solution guess. The segment signature, contained in
% the modes, events, and resets variables uniquely defines the multisegment
% boundary-value problem. Continuation is restricted to two-segment
% periodic orbits without impacts and with stick throughout.

p0 = [0.59; 0.04; 3.17; 0.92; 0.80];
modes  = {'stick' 'stick'};
events = {'phase' 'minsep'};
resets = {'phase' 'turn'};
f  = @(t, x) stickslip(x, p0, modes{1});
[t1, x1] = ode45(f, [0  1.5], [0.5; 0; 0]);
f  = @(t, x) stickslip(x, p0, modes{2});
[t2, x2] = ode45(f, [0  1.5], [0.25; 0; -pi/2]);
t0 = {t1 t2};
x0 = {x1 x2};

% Turn off bifurcation detection

prob = coco_prob();
prob = coco_set(prob, 'hspo', 'bifus', false);
prob = ode_isol2hspo(prob, '', ...
  {@stickslip, @stickslip_events, @stickslip_resets}, ...
  modes, events, resets, t0, x0, {'V' 'c' 'n' 'w' 'e'}, p0);

fprintf('\n Run=''%s'': Continue family of two-segment periodic orbits with sustained stick.\n', ...
  'stickslip1');

coco(prob, 'stickslip1', [], 1, 'V', [0.5 0.7]);

%% Locate grazing periodic orbit

% The continuation problem structure encoded below includes an additional
% monitor function that evaluates to the first component of the initial end
% point on the first trajectory segment, with corresponding inactive
% continuation parameter 'pos', whose initial value is assigned to equal
% 0.5. Its dimensional deficit equals -1. A unique solution to a closed
% continuation problem is obtained by releasing 'V', while keeping 'pos'
% fixed. The corresponding periodic orbit grazes the impact surface.

[data, uidx] = coco_get_func_data(prob, 'hspo.orb.bvp', 'data', 'uidx');
prob = coco_add_pars(prob, 'graze', uidx(data.bvp_bc.x0_idx(1)), 'pos');
prob = coco_set_parival(prob, 'pos', 0.5);

fprintf(...
  '\n Run=''%s'': Locate a grazing periodic orbit.\n', ...
  'graze_run');

coco(prob, 'graze_run', [], 0, {'V' 'pos'});

%% Continue family of stick-slip periodic orbits with impacts

% The continuation problem structure encoded below includes five monitor
% functions that evaluate to the problem parameters, and the corresponding
% inactive continuation parameters 'V', 'c', 'n', 'w', and 'e'. A
% one-dimensional family of periodic, impacting stick-slip oscillations
% results by releasing 'V' and allowing it to vary during continuation.

prob = coco_prob();
prob = coco_set(prob, 'hspo', 'bifus', false);

% Construct an initial solution guess by reading the grazing periodic orbit
% from the previous run and appending a zero-length segment. The associated
% changes to segment signature impose boundary conditions corresponding to
% transitions between stick to slip due to impacts and from slip back to
% stick when the slipping velocity equals 0.

modes  = {'stick' 'stick' 'slip'};
events = {'phase' 'collision' 'rest'};
resets = {'phase' 'bounce' 'stick'};
sol = hspo_read_solution('', 'graze_run', 1);
t0 = {sol.tbp{1}, sol.tbp{2}, 0};                   % Third segment has zero length
x0 = {sol.xbp{1}, sol.xbp{2}, [0 sol.xbp{1}(1,:)]}; % And consists of a single point
p0 = sol.p;

prob = ode_isol2hspo(prob, '', ...
  {@stickslip, @stickslip_events, @stickslip_resets}, ...
  modes, events, resets, t0, x0, {'V' 'c' 'n' 'w' 'e'}, p0);

fprintf('\n Run=''%s'': Continue family of three-segment periodic orbits with impacts and stick-slip transitions.\n', ...
  'stickslip2');

coco(prob, 'stickslip2', [], 1, 'V', [0.5 0.7]);

%% Graphical representation of results

% Plot data: panels (a)-(c)
labs = [5, 1, 3];
seg1.sol.RO = {'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.5 0.5 0.5], 'Marker', '.', 'MarkerSize', 12};
seg2.sol.RO = {'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black'};
for i=1:3
  figure(i); clf; hold on; grid on; box on; axis([0 0.8 -0.3 0.3])
  
  plot([0.5 0.5], [-0.3 0.3], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', 'black');
  coco_plot_sol(seg1', 'stickslip1', labs(i), 'hspo.orb.bvp.seg1', 'x', 'x')
  coco_plot_sol(seg2', 'stickslip1', labs(i), 'hspo.orb.bvp.seg2', 'x', 'x')
  
  hold off
end

% Plot data: panels (a)-(c)
labs = [4, 1, 2];
seg1.sol.RO = {'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.5 0.5 0.5], 'Marker', '.', 'MarkerSize', 12};
seg2.sol.RO = {'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black'};
seg3.sol.RO = {'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black'};
for i=1:3
  figure(i+3); clf; hold on; grid on; box on; axis([0 0.8 -0.4 0.4])
  
  plot([0.5 0.5], [-0.4 0.4], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', 'black')
  coco_plot_sol(seg1', 'stickslip2', labs(i), 'hspo.orb.bvp.seg1', 'x', 'x')
  coco_plot_sol(seg2', 'stickslip2', labs(i), 'hspo.orb.bvp.seg2', 'x', 'x')
  coco_plot_sol(seg2', 'stickslip2', labs(i), 'hspo.orb.bvp.seg3', 'x', 'x')
  
  hold off
end
