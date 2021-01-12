%% Bang-bang excitation of the Duffing oscillator from Examples 17.3 in Recipes for Continuation
%
% We investigate families of periodic orbits of the nonlinear Duffing
% oscillator under bang-bang excitation. The analysis is most conveniently
% expressed in terms of multi-segment boundary-value problems with
% individual segments characterized by one of the two values of the
% excitation.

% Panel (a) shows a primary branch of two-segment periodic orbits along
% which black indicates asymptotically stable orbits and grey represents
% unstable orbits. Several secondary branches of two-segment periodic
% orbits emanating from branch points along the primary branch are included
% in the figure, together with branches of four-segment periodic orbits
% emanating from period-doubling bifurcations along the secondary branches.
% The figure also includes saddle-node and period-doubling bifurcation
% curves that intersect the primary and secondary branches, respectively.
% Panels (b), (c), and (d) shows sample two-segment periodic orbits along
% the primary branch.
 
%% Encoding

% The continuation problem structure encoded below includes five monitor
% functions that evaluate to the problem parameters, and the corresponding
% inactive continuation parameters 'la', 'al', 'eps', 'A', and 'om'. Its
% dimensional deficit equals 0. A one-dimensional family of periodic orbits
% is obtained by releasing 'om' and allowing it vary during continuation.
% The atlas algorithm detects fold points and branch points and identifies
% these by the event types 'FP' and 'BP', respectively. By default, special
% points associated with saddle-node bifurcations and period-doubling
% bifurcations are detected by the 'hspo' toolbox and identified by the
% event types 'SN' and 'PD'. Notably, 'SN' events are also sometimes
% triggered at 'BP' points.

% Construct initial solution guess
w0 = 1.1;
p0 = [0.2; 1; 1; 26; w0];
x0 = [0; 0; 0];
for i=1:10 % Let transients die out
  [~, y0] = ode45(@(t,x) duff(x,p0,'neg'), [0 pi/w0], x0);
  x0 = y0(end,:)';
  [~, y1] = ode45(@(t,x) duff(x,p0,'pos'), [0 pi/w0], x0);
  x0 = y1(end,:)';
end

x0 = [x0(1:2);0];
[t1, x1] = ode45(@(t,x) duff(x,p0,'neg'), [0 pi/w0], x0); % First segment of approximate periodic orbit
x0 =  [x1(end,1:2) 0]';
[t2, x2] = ode45(@(t,x) duff(x,p0,'pos'), [0 pi/w0], x0); % Second segment of approximate periodic orbit

modes  = {'neg' 'pos'};
events = {'phase' 'phase'};
resets = {'phase' 'phase'};
t0     = {t1 t2};
x0     = {x1 x2};

% Use anonymous, in-line definitions of the event and reset functions and
% their Jacobians.

duff_events    = @(x,p,~) pi/p(5)-x(3);
duff_events_dx = @(x,p,~) [0 0 -1];
duff_events_dp = @(x,p,~) [0 0 0 0 -pi/p(5)^2];
duff_resets    = @(x,p,~) [x(1);x(2);0];
duff_resets_dx = @(x,p,~) [1 0 0; 0 1 0; 0 0 0];
duff_resets_dp = @(x,p,~) zeros(3,5);

% Set a high value for the NTST discretization parameter for all 'coll'
% instances.

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 200);
prob = ode_isol2hspo(prob, '', {@duff, duff_events, duff_resets}, ...
  {@duff_DFDX, duff_events_dx, duff_resets_dx}, ...
  {@duff_DFDP, duff_events_dp, duff_resets_dp}, ...
  modes, events, resets, t0, x0, {'la' 'al' 'eps' 'A' 'om'}, p0);

% The @duff_add_IP slot function responds to the 'bddat' signal and stores
% the components of the initial end point on the first trajectory segment
% to the cell array returned by the coco entry-point function and stored to
% disk.

prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 20, 'PtMX', 200, 'NAdapt', 1);

fprintf('\n Run=''%s'': Continue primary family of two-segment periodic orbits.\n', ...
  'run1');

coco(prob, 'run1', [], 1, 'om', [0.8 1.5]);

%% Continue along secondary branches of two-segment periodic orbits

% The encoding below is identical to the one above.

bd1  = coco_bd_read('run1');
labs = coco_bd_labs(bd1, 'BP');
for lab=labs
  prob = coco_prob();
  prob = ode_BP2hspo(prob, '', 'run1', lab);
  prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
  prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 10, 'PtMX', 50, 'NAdapt', 5);
  
  fprintf(...
  '\n Run=''%s'': Continue secondary family of two-segment periodic orbits from point %d in run ''%s''.\n', ...
  sprintf('run2_%d',lab), lab, 'run1');

  coco(prob, sprintf('run2_%d',lab), [], 1, 'om', [0.5 1.5]);
end

%% Continuation of four-segment periodic orbits

% The encoding below uses solution data stored with a period-doubling
% bifurcation located along a secondary branch to restart continuation of a
% four-segment periodic orbit. The dimensional deficit is again 0,
% requiring the release of a single continuation parameter in order to
% obtain a one-dimensional solution manifold.

run = sprintf('run2_%d',labs(1));
bd2 = coco_bd_read(run);
labs = coco_bd_labs(bd2, 'PD');
for lab=labs([1 end])
  prob = coco_prob();
  prob = ode_PD2hspo(prob, '', run, lab);
  prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
  prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 10, 'PtMX', 50, 'NAdapt', 5);
  
  fprintf(...
    '\n Run=''%s'': Continue family of period-doubled four-segment periodic orbits from point %d in run ''%s''.\n', ...
    sprintf('run3_%d', lab), lab, run);

  coco(prob, sprintf('run3_%d', lab), [], 1, 'om', [0.5 1.5]);
end

%% Continuation of period-doubling bifurcations

% The encoding below imposes additional constraints on the two-segment
% boundary-value problem corresponding to a period-doubling bifurcation,
% resulting in a dimensional deficit of -1. A one-dimensional solution
% manifold results by releasing 'om' and 'A' and allow these to vary during
% continuation.

bd1  = coco_bd_read('run1');
labs = coco_bd_labs(bd1, 'BP');
run  = sprintf('run2_%d',labs(1));
bd2  = coco_bd_read(run);
labs = coco_bd_labs(bd2, 'PD');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 20, 'PtMX', 100, 'NAdapt', 1);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
prob = ode_PD2PD(prob, '', run, labs(1));

fprintf(...
  '\n Run=''%s'': Continue family of period-doubling bifurcations from point %d in run ''%s''.\n', ...
  'run4', labs(1), run);

coco(prob, 'run4', [], 1, {'om' 'A'}, [0.5 1.5]);

%% Continuation of saddle-node bifurcations

% The encoding below imposes additional constraints on the two-segment
% boundary-value problem corresponding to a saddle-node bifurcation,
% resulting in a dimensional deficit of -1. A one-dimensional solution
% manifold results by releasing 'om' and 'A' and allow these to vary during
% continuation.

bd1  = coco_bd_read('run1');
labs = coco_bd_labs(bd1, 'SN');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 20, 'PtMX', 200, 'NAdapt', 5);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
prob = ode_SN2SN(prob, '', 'run1', labs(1));

fprintf(...
  '\n Run=''%s'': Continue family of saddle-node bifurcations from point %d in run ''%s''.\n', ...
  'run5', labs(1), 'run1');

coco(prob, 'run5', [], 1, {'om' 'A'}, [0.2 1.5]);

%% Graphical representation of results

% Panel (a)
bd1 = coco_bd_read('run1');

figure(1); clf; hold on; grid on; box on; view(160,25)

thm = struct('special', {{'EP', 'FP', 'BP', 'PD'}});
coco_plot_bd(thm, 'run1', 'om', 'A', 'X0')

labs = coco_bd_labs(bd1, 'BP');

for lab=labs
  coco_plot_bd(thm, sprintf( 'run2_%d',lab), 'om', 'A', 'X0')
end

run = sprintf('run2_%d',labs(1));
bd2 = coco_bd_read(run);
labs = coco_bd_labs(bd2, 'PD');

for lab=labs([1 end])
  coco_plot_bd(thm, sprintf('run3_%d',lab), 'om', 'A', 'X0')
end

thm = struct('special', {{'EP', 'FP'}}, ...
  'xlab', '\omega', 'zlab', 'x_1(0)');
coco_plot_bd(thm, 'run4', 'om', 'A', 'X0')
coco_plot_bd(thm, 'run5', 'om', 'A', 'X0')

hold off

% Plot data: panels (b)-(d)
labs  = [9 12 16];

seg1.sol.RO = {'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12};
seg2.sol.RO = {'LineStyle', '-', 'LineWidth', 2, 'Color', ...
  [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12};

for i=1:3
  figure(i+1); clf; hold on; box on; grid on
  coco_plot_sol(seg1, 'run1', labs(i), 'hspo.orb.bvp.seg1', 'x', 'x')
  coco_plot_sol(seg2, 'run1', labs(i), 'hspo.orb.bvp.seg2', 'x', 'x')
  axis([-inf inf -inf inf]); hold off
end
