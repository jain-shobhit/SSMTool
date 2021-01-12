%% One-dimensional manifold continuation using the atlas_kd atlas algorithm

coco_func_data.pointers('set', []);

prob = coco_prob();
prob = coco_set(prob, 'cont', 'atlas', 'kd', 'PtMX', 250);
prob = coco_add_func(prob, 'ellipse', @ellipse, [], ...
  'zero', 'u0', [1.1; 0]);

%% Cartesian coordinates, no computational boundary (stops)

pprob = coco_add_pars(prob, 'cartpars', 1:2, {'x', 'y'}, 'active');
coco(pprob, 'run1', [], 1, {'x' 'y'})

atlas = coco_bd_read('run1', 'atlas');
figure(1); clf
plot_atlas_kd(atlas.charts, 2)
axis equal; hold on
t=0:0.01:2*pi;
plot(1+cos(t)/10,sin(t),'b-.')
hold off

%% Cartesian coordinates, with computational domain boundary (stops)
coco(pprob, 'run2', [], 1, {'x' 'y'}, {[] [0 1]})

atlas = coco_bd_read('run2', 'atlas');
figure(2); clf
plot_atlas_kd(atlas.charts, 2)
axis tight; axis equal

%% Cartesian coordinates and active angle parameter (doesn't stop)

pprob = coco_add_pars(prob, 'cartpars', [1 2], {'x' 'y'});
pprob = coco_add_func(pprob, 'arclength', @arclength, [], ...
  'zero', 'uidx', [1 2], 'u0', 0);
pprob = coco_add_pars(pprob, 'phipar', 3, 'phi');
coco(pprob, 'run3', [], 1, {'x' 'y' 'phi'})

atlas = coco_bd_read('run3', 'atlas');
figure(3); clf
plot_atlas_kd(atlas.charts, 2)
axis tight; axis equal

%% One Cartesian coordinate and active angle parameter (doesn't stop)

pprob = coco_add_pars(prob, 'cartpar', 1, 'x');
pprob = coco_add_func(pprob, 'arclength', @arclength, [], ...
  'zero', 'uidx', [1 2], 'u0', 0);
pprob = coco_add_pars(pprob, 'phipar', 3, 'phi');
coco(pprob, 'run4', [], 1, {'x' 'phi'})

atlas = coco_bd_read('run4', 'atlas');
figure(4); clf
plot_atlas_kd(atlas.charts, 2)
axis normal

pprob = coco_add_pars(prob, 'cartpar', 2, 'y');
pprob = coco_add_func(pprob, 'arclength', @arclength, [], ...
  'zero', 'uidx', [1 2], 'u0', 0);
pprob = coco_add_pars(pprob, 'phipar', 3, 'phi');
coco(pprob, 'run5', [], 1, {'y' 'phi'})

atlas = coco_bd_read('run5', 'atlas');
figure(5); clf
plot_atlas_kd(atlas.charts, 2)
axis normal

%% Graphical representation of stored solutions

phi=0:2*pi/100:2*pi;
% Terminate near singular projection onto x-axis
figure(6); clf; hold on
plot(1+cos(phi)/10,sin(phi),'b--')
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
thm.special = {'MX'};
thm.MX = {'d', 'MarkerSize', 8, 'MarkerFaceColor', 'k'};
coco_plot_bd(thm, 'run1')
axis equal; grid on; hold off

% Projection onto angular variable
figure(7); clf; hold on
plot(1+cos(phi)/10,sin(phi),'b--')
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'run3', 'x', {'x' 'phi'}, @(x,phi) (x-1).*tan(phi))
axis equal; grid on; hold off

% Projection onto (x,y)-plane
figure(8); clf; hold on
plot(1+cos(phi)/10,sin(phi),'b--')
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'run4', 'x', {'x' 'phi'}, @(x,phi) (x-1).*tan(phi))
axis equal; grid on; hold off

% Projection onto (x,y)-plane
figure(9); clf; hold on
plot(1+cos(phi)/10,sin(phi),'b--')
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'run5', {'y' 'phi'}, @(y,phi) 1+y.*cot(phi), 'y')
axis equal; grid on; hold off
