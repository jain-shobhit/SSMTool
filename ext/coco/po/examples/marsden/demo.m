%% A Hopf bifurcation from Example 8.3 in Recipes for Continuation
%
% We compute a family of periodic orbits emanating from a Hopf bifurcation
% point of the dynamical system
%
% x1' = p1*x1 + x2 + p2*x1^2
% x2' = -x1 + p1*x2 + x2*x3
% x3' = (p1^2 - 1)*x2 - x1 - x3 + x1^2

% We obtain an initial solution guess from normal form analysis, shown in
% gray in panel (a). The initial correction step coverges to orbit 1. Panel
% (b) shows the family of periodic orbits of increasing amplitudes that
% seem to approach a homoclinic orbit, indicated by the corner that
% develops in the top left part of the plot and allocates many mesh points
% due to slow dynamics.
%
% Panels (c) and (d) show the last periodic orbit obtained in the
% continuation run in state space and as a time history. The time history
% shows a phase of slow dynamics, indicating existence of a nearby
% equilibrium point, followed by a fast excursion. 
%
% Starting with the last periodic orbit, we insert a long segment
% of constant dynamics and rescale the period such that the shape of the
% orbit in phase space should be unchanged if there exists a nearby
% homoclinic orbit. The extended time profile after the initial correction
% step is shown in panel (e). We clearly observe an elongated phase of
% near-constant dynamics. We overlay this new solution (black dot) on top
% of the previous orbit (gray circle) in panel (f). The phase plots,
% including the distribution of mesh points, are virtually identical, which
% supports the assumption that a nearby homoclinic orbit exists.
%
% Continuation of the periodic orbit with high period, while keeping the
% period constant, resulting in an approximation to a homoclinic
% bifurcation curve (g). Each point on this curve corresponds to a terminal
% point along a family of periodic orbits emenating from a Hopf bifurcation
% under variations in p_2. Panel (h) shows selected members of the family
% of high-period orbits.

%% Encoding

% The continuation problem encoded below includes three monitor functions
% that evaluate to the problem parameters and the interval length,
% respectively, and the corresponding inactive continuation parameters 'p1'
% and 'p2' and active continuation parameter 'po.period'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'p1' is released and allowed to vary during continuation.

prob = coco_prob();
prob = ode_isol2ep(prob, '', @marsden, [0; 0; 0], {'p1', 'p2'}, [-1; 6]);

fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep_run');

bd1  = coco(prob, 'ep_run', [], 1, 'p1', [-1 1]);

HBlab = coco_bd_labs(bd1, 'HB');

prob = coco_prob();
prob = ode_HB2po(prob, '', 'ep_run', HBlab);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', 50);

fprintf(...
  '\n Run=''%s'': Continue family of periodic orbits from point %d in run ''%s''.\n', ...
  'po1', HBlab, 'ep_run');

bd2  = coco(prob, 'po1', [], 1, {'p1' 'po.period'}, [-1 1]);

%% Locate a high-period periodic orbit that approximates a homoclinic orbit

% The continuation problem structure encoded below is identical to that
% above, but an initial solution guess is obtained by extending the
% duration spent near an equilibrium for a periodic orbit found in the
% previous run. The call to coco_xchg_pars constrains the interval length,
% while releasing the second problem parameter.

[sol, data] = coll_read_solution('po.orb', 'po1', 8); % Periodic orbit with T=1.8652e+01
f = marsden(sol.xbp', repmat(sol.p, [1 size(sol.xbp, 1)])); % Evaluate vector field at basepoints
[~, idx] = min(sqrt(sum(f.*f, 1))); % Find basepoint closest to equilibrium

scale = 25;
T  = sol.T;
t0 = [sol.tbp(1:idx,1) ; T*(scale-1)+sol.tbp(idx+1:end,1)]; % Crank up period by factor scale
x0 = sol.xbp;
p0 = sol.p;

% Initialize continuation problem structure with the same number of
% intervals as in previous run.

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', data.coll.NTST);
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_isol2po(prob, '', @marsden, t0, x0, {'p1' 'p2'}, p0);
prob = coco_set(prob, 'cont', 'NAdapt', 10);
prob = coco_xchg_pars(prob, 'p2', 'po.period');

fprintf(...
  '\n Run=''%s'': Find reconstructed high-period periodic orbit approximating a homoclinic connection.\n', ...
  'po2');

coco(prob, 'po2', [], 0, {'p1' 'po.orb.coll.err_TF' 'po.period'});

% The family of high-period periodic orbits of constant period
% approximates a family of homoclinic connections to a saddle equilibrium.

prob = coco_prob();
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_po2po(prob, '', 'po2', 2);
prob = coco_xchg_pars(prob, 'p2', 'po.period');
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', 50);

fprintf(...
  '\n Run=''%s'': Continue family of high-period periodic orbits from point %d in run ''%s''.\n', ...
  'po3', 2, 'po2');

coco(prob, 'po3', [], 1, {'p1' 'p2' 'po.period'}, [-1 1]);

%% Graphical representation of stored solutions

% Plot data: panel (a)

figure(1); clf; hold on; grid on; box on; axis(0.001*[-0.8 0.8 -0.8 0.8])
t0 = (0:2*pi/100:2*pi)';
x0 = 0.001*(cos(t0)*[1 0 -1]-sin(t0)*[0 1 0])/sqrt(2);
plot(x0(:,2), x0(:,3), 'ro')
coco_plot_sol('po1', 1, '', 'x', 2, 'x', 3)
hold off

% Plot data: panel (b)

figure(2); clf; 
coco_plot_sol('po1', 1:8, '', 'x', 2, 'x', 3)
box on; grid on; axis([-0.5 0.1 -0.1 0.7])

% Plot data: panel (c)

figure(3); clf; 
coco_plot_sol('po1', 8, '', 'x', 2, 'x', 3)
grid on; box on; axis([-0.5 0.1 -0.1 0.7])

% Plot data: panel (d)

figure(4); clf; hold on; grid on; box on; axis([-0.5 20 -0.1 0.8])

sol = po_read_solution('po1', 8); % Extract solution
f = marsden(sol.xbp', repmat(sol.p, [1 size(sol.xbp, 1)]));
[~, idx] = min(sqrt(sum(f.*f, 1))); % Find basepoint closest to equilibrium
tbp = [sol.tbp(idx:end)-sol.tbp(idx); ...
  sol.tbp(1:idx)+sol.T-sol.tbp(idx)];
xbp = sol.xbp([idx:end 1:idx],3);
plot(tbp, xbp, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
xlabel 't'
ylabel 'x_3'

hold off

% Plot data: panel (e)

figure(5); clf; hold on; grid on; box on; axis([-10 500 -0.1 0.8])

sol1 = po_read_solution('po2', 1); % Extract solution
f = marsden(sol1.xbp', repmat(sol1.p, [1 size(sol1.xbp, 1)]));
[~, idx] = min(sqrt(sum(f.*f, 1))); % Find basepoint closest to equilibrium
tbp = [sol1.tbp(idx:end)-sol1.tbp(idx); ...
  sol1.tbp(1:idx)+sol1.T-sol1.tbp(idx)];
xbp = sol1.xbp([idx:end 1:idx],3);
plot(tbp, xbp, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
xlabel 't'
ylabel 'x_3'

hold off

% Plot data: panel (f)

figure(6); clf; hold on; grid on; box on; axis([-0.5 0.1 -0.1 0.7])
coco_plot_sol('po1', 8, '', 'x', 2, 'x', 3)
coco_plot_sol('po2', 1, '', 'x', 2, 'x', 3)
hold off

% Plot data: panel (g)

figure(7); clf
thm = struct('ustab', '', 'xlab', 'p_1', 'ylab', 'p_2');
thm.lspec = {'k', 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 8};
coco_plot_bd(thm, 'po3', 'p1', 'p2')
grid on; box on; axis([-0.08 0.0 1 21])

% Plot data: panel (h)

figure(8); clf 
coco_plot_sol('po3', [1 4 6 9 12], '', 'x', 2, 'x', 3)
grid on; box on; axis([-0.7 0.4 -0.3 1.5])
