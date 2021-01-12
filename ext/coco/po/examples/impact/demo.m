%% An impact oscillator from Example 15.2 in Recipes for Continuation
%
% We investigate families of impacting periodic orbits of a forced, linear
% mechanical oscillator whose motion is constrained by a mechanical stop.
% Contact is modeled as instantaneous (sustained contact is not considered)
% and associated with a velocity jump governed by a Newtonian coefficient
% of restitution.

% Panel (a) shows the approximate L2 norm along a family of impacting
% periodic orbits under variations in the excitation amplitude A. Panels
% (b), (c), and (d) show three solutions to the corresponding two-segment
% boundary-value problem, but the solution in (d) is not physically
% realizable, since it includes excursion beyond the mechanical constraint.
% Panel (e) shows a two-parameter representation of a family of two-segment
% grazing periodic orbits. Corresponding sample orbits are shown in panels
% (f), (g), and (h), all of which are physically realizable.

%% Initial encoding

% The continuation problem structure encoded below includes six monitor
% functions that evaluate to the problem parameters and a single monitor
% function that evaluates to the second component (the physical velocity)
% of the end point at t=0 of the second trajectory segment, and the
% corresponding inactive continuation parameters 'k', 'c', 'A', 'w', 'd',
% and 'e', and active continuation parameter 'graze'. Its dimensional
% deficit equals 0. A one-dimensional family of impacting periodic orbits
% is obtained by releasing 'A' and allowing this to vary during
% continuation. Events associated with zero-crossings of the value of
% 'graze' (grazing periodic orbits) are detected and identified with the
% 'GR' event type.

% Construct initial solution guess. The segment signature, contained in
% the modes, events, and resets variables uniquely defines the multisegment
% boundary-value problem. Continuation is restricted to two-segment
% periodic orbits with impact surface at x1=d.

p0     = [1; 0.1; 1; 1; 1; 0.8];
modes  = {'free' 'free'};
events = {'impact' 'phase'};
resets = {'bounce' 'phase'};
f       = @(t, x) impact(x, p0, 'free');
[t1, x1] = ode45(f, [0 3.2], [-0.98; -0.29; -pi]);
[t2, x2] = ode45(f, [0 3.1], [1; -1.36; 0.076]);
t0 = {t1  t2};
x0 = {x1  x2};
hspo_args = {{@impact, @impact_events, @impact_resets}, ...
  modes, events, resets, t0, x0, {'k' 'c' 'A' 'w' 'd' 'e'}, p0};

% Turn off bifurcation detection

prob = coco_prob();
prob = coco_set(prob, 'hspo', 'bifus', false);
prob = ode_isol2hspo(prob, '', hspo_args{:});

% Extract function data and dependency index set for second trajectory
% segment and use this to extract the second component of the initial end
% point on this segment. The active continuation parameter 'graze' tracks
% the value of this component, corresponding to the velocity immediately
% after impact.

[data, uidx] = coco_get_func_data(prob, 'hspo.orb.bvp.seg2.coll', ...
  'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'grazing', uidx(maps.x0_idx(2)), ...
  'graze', 'active');
prob = coco_add_event(prob, 'GR', 'graze', 0);
prob = coco_set(prob, 'cont', 'PtMX', 100, 'NAdapt', 5);

fprintf('\n Run=''%s'': Continue family of two-segment impacting periodic orbits.\n', ...
  'impact1');

bd1 = coco(prob, 'impact1', [], {'A' 'graze'}, [0.01 1]);


%% Continuation of grazing periodic orbits

% The continuation problem structure encoded below is identical to that
% above, but constructed from stored data associated with the grazing
% periodic orbit, and with a reassignment of 'graze' to the inactive
% continuation parameters and 'A' to the active continuation parameters. A
% one-dimensional family of grazing periodic orbits is obtained by
% releasing 'w' and allowing it to vary during continuation.

labgr = coco_bd_labs(bd1, 'GR');
prob = coco_prob();
prob = coco_set(prob, 'hspo', 'bifus', false);
prob = ode_hspo2hspo(prob, '', 'impact1', labgr);
[data, uidx] = coco_get_func_data(prob, 'hspo.orb.bvp.seg2.coll', ...
  'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'grazing', uidx(maps.x0_idx(2)), ...
  'graze', 'active');

% Switch activation status between 'graze' (initially active) and 'A'
% (initially inactive).

prob = coco_xchg_pars(prob, 'graze', 'A');
prob = coco_set(prob, 'cont', 'PtMX', 100, 'NAdapt', 5);

fprintf(...
  '\n Run=''%s'': Continue grazing periodic orbits from point %d in run ''%s''.\n', ...
  'impact2', labgr, 'impact1');

coco(prob, 'impact2', [], {'w' 'A' 'graze'}, {[] [0 1]});

%% Graphical representation of results

figure(1); clf;
thm = struct('ustab', '', 'special', {{'EP', 'GR'}});
thm.EP = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'cyan', 'MarkerFaceColor', 'white'};
thm.GR = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white'};
thm.ylab = 'L_2[0,T]';
coco_plot_bd(thm, 'impact1', 'A', ...
  {'||hspo.orb.bvp.seg1.x||_{L_2[0,T]}', ...
  '||hspo.orb.bvp.seg2.x||_{L_2[0,T]}'}, ...
  @(x,y) sqrt(x.^2+y.^2))
grid on; box on; axis([0.05 1.05 5.1 6])

labs = [1, 6, 11];
seg1.sol.RO = {'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12};
seg2.sol.RO = {'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12};
for i=1:3
  figure(i+1); clf; hold on; grid on; box on; axis([-1.5 2 -2 2])
  
  coco_plot_sol(seg1, 'impact1', labs(i), 'hspo.orb.bvp.seg1', 'x', 'x')
  coco_plot_sol(seg2, 'impact1', labs(i), 'hspo.orb.bvp.seg2', 'x', 'x')
  plot([1 1], [-2 2], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', 'black')
  xlabel('x_1'); ylabel('x_2')
  
  hold off
end

figure(5); clf; hold on
thm = struct('ustab', '', 'special', {{'EP'}}, 'ylab', '\omega');
thm.EP = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'cyan', 'MarkerFaceColor', 'white'};
coco_plot_bd(thm, 'impact2', 'A', 'w')
grid on; box on; axis([0 1 0 1.5])

labs = [16, 1, 11];

for i=1:3
  figure(i+5); clf; hold on
  
  coco_plot_sol(seg1, 'impact2', labs(i), 'hspo.orb.bvp.seg1', 'x', 'x')
  coco_plot_sol(seg2, 'impact2', labs(i), 'hspo.orb.bvp.seg2', 'x', 'x')
  plot([1 1], [-2 2], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', 'black')
  xlabel('x_1'); ylabel('x_2');
  
  hold off; grid on; box on; axis([-1.1 1.1 -2 2])
end
