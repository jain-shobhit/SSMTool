addpath ../atlas2d_v6 % Add atlas2d_v6 atlas algorithm to search path

prob = coco_prob();
prob = coco_set(prob, 'cont', 'linsolve', 'recipes');  % linear solver
prob = coco_set(prob, 'cont', 'corrector', 'recipes'); % nonlinear solver
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'all', 'CleanData', true);

%% Tests of performance for different choices of projected geometry

pprob = coco_add_func(prob, 'torus', @torus, [], ...
  'zero', 'u0', [0; 2; 1]);
pprob = coco_add_pars(pprob, 'cartpar', [1 2 3], {'x' 'y' 'z'});

pprob = coco_set(pprob, 'cont', 'R', 0.3, 'PtMX', 500, ...
  'indcs', [1 2 3], 'almax', 30);
coco(pprob, 'torus1', [], 2, {'x' 'y' 'z'}, [0 2])

pprob = coco_add_func(pprob, 'angles', @angles, [], ...
  'zero', 'uidx', [1 2 3], 'u0', [pi/2; pi/2]);
pprob = coco_add_pars(pprob, 'phipar', [4 5], {'phi' 'theta'});

pprob = coco_set(pprob, 'cont', 'R', 0.4, 'PtMX', 500, 'indcs', [4 5]);
coco(pprob, 'torus2', [], 2, {'x' 'y' 'z' 'phi' 'theta'}, [0 2])

pprob = coco_set(pprob, 'cont', 'R', 0.4, 'PtMX', 500, 'indcs', [1 2 3]);
coco(pprob, 'torus3', [], 2, {'x' 'y' 'z' 'phi' 'theta'}, [0 2])

%% Graphical representation of stored solutions

n = 20;
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi;
cosphi = cos(phi);
sintheta = sin(theta);

X = (2+0.99.*cosphi)*cos(theta);
Y = (2+0.99.*cosphi)*sin(theta);
Z = 0.99.*sin(phi)*ones(1,n+1);

% Projection onto (x,y,z)-space without angles.
figure(1); clf; hold on
surf(X,Y,Z, 'EdgeColor', 0.7*[1 1 1], 'FaceColor', 0.8*[1 1 1]);
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'torus1', 'x', 'y', 'z')
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

figure(2); clf; hold on
atlas = coco_bd_read('torus1', 'atlas');
plot_charts(atlas.charts, 1, 2, 3)
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

figure(3); clf; hold on
plot_trisurf(atlas.charts, 1, 2, 3)
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

% Projection onto angular variables
figure(4); clf; hold on
surf(X,Y,Z, 'EdgeColor', 0.7*[1 1 1], 'FaceColor', 0.8*[1 1 1]);
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'torus2', 'x', 'y', 'z')
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

figure(5); clf; hold on
atlas = coco_bd_read('torus2', 'atlas');
plot_charts(atlas.charts, 1, 2, 3)
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

figure(6); clf; hold on
plot_trisurf(atlas.charts, 1, 2, 3)
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

% Projection onto (x,y,z)-space with angles
figure(7); clf; hold on
surf(X,Y,Z, 'EdgeColor', 0.7*[1 1 1], 'FaceColor', 0.8*[1 1 1]);
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'torus3', 'x', 'y', 'z')
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

figure(8); clf; hold on
atlas = coco_bd_read('torus3', 'atlas');
plot_charts(atlas.charts, 1, 2, 3)
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

figure(9); clf; hold on
plot_trisurf(atlas.charts, 1, 2, 3)
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

rmpath ../atlas2d_v6 % Remove atlas2d_v5 atlas algorithm from search path
