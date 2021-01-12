addpath ../atlas1d_v5 % Add atlas1d_v5 atlas algorithm to search path

prob = coco_prob();
prob = coco_set(prob, 'cont', 'linsolve', 'recipes');  % linear solver
prob = coco_set(prob, 'cont', 'corrector', 'recipes'); % nonlinear solver
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'all', 'CleanData', true);

%% Tests of performance for different choices of projected geometry

pprob = coco_add_func(prob, 'circle', @circle, [], ...
  'zero', 'u0', [1.5; 1]);
pprob = coco_add_pars(pprob, 'cartpar', [1 2], {'x' 'y'});

pprob = coco_set(pprob, 'cont', 'R', 0.1, 'PtMX', -30, 'indcs', 1);
coco(pprob, 'run1', [], 1, {'x' 'y'}, [1.5 2])

pprob = coco_set(pprob, 'cont', 'R', 0.1, 'PtMX', 30, 'indcs', 2);
coco(pprob, 'run2', [], 1, {'x' 'y'}, [1.5 2])

pprob = coco_add_func(pprob, 'angle', @angle, [], ...
  'zero', 'uidx', [1 2], 'u0', atan2(1,1.5));
pprob = coco_add_pars(pprob, 'phipar', 3, 'phi');

pprob = coco_set(pprob, 'cont', 'R', 0.1, 'PtMX', -30, 'indcs', 3);
coco(pprob, 'run3', [], 1, {'x' 'y' 'phi'}, [1.5 2])

pprob = coco_set(pprob, 'cont', 'R', 0.1, 'PtMX', -30, 'indcs', [1 2]);
coco(pprob, 'run4', [], 1, {'x' 'y' 'phi'}, [1.5 2])

rmpath ../atlas1d_v5 % Remove atlas1d_v4 atlas algorithm from search path

%% Graphical representation of stored solutions

phi=0:2*pi/100:2*pi;

% Gap near singular projection onto x-axis
figure(1); clf; hold on
plot(1+cos(phi),sin(phi),'b--')
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
thm.special = {'MX'};
thm.MX = {'d', 'MarkerSize', 8, 'MarkerFaceColor', 'k'};
coco_plot_bd(thm, 'run1')
coco_plot_bd('run1')
axis equal; grid on; hold off

% Immediate termination due to gap
figure(2); clf; hold on
plot(1+cos(phi),sin(phi),'b--')
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
thm.special = {'MX'};
thm.MX = {'d', 'MarkerSize', 8, 'MarkerFaceColor', 'k'};
coco_plot_bd(thm, 'run2')
coco_plot_bd('run2')
axis equal; grid on; hold off

% Projection onto angular variable
figure(3); clf; hold on
plot(1+cos(phi),sin(phi),'b--')
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'run3')
coco_plot_bd('run3')
axis equal; grid on; hold off

% Projection onto (x,y)-plane
figure(4); clf; hold on
plot(1+cos(phi),sin(phi),'b--')
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'run4')
coco_plot_bd('run4')
axis equal; grid on; hold off
