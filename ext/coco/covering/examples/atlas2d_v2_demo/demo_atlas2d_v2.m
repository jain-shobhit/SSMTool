addpath ../atlas2d_v2 % Add atlas2d_v2 atlas algorithm to search path

prob = coco_prob();
prob = coco_set(prob, 'cont', 'linsolve', 'recipes');  % linear solver
prob = coco_set(prob, 'cont', 'corrector', 'recipes'); % nonlinear solver
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'all', 'CleanData', true);

%% Tests of performance for different choices of projected geometry

pprob = coco_add_func(prob, 'torus', @torus, [], ...
  'zero', 'u0', [1; -1; -1]);
pprob = coco_add_pars(pprob, 'cartpar', [1 2 3], {'x' 'y' 'z'});

pprob = coco_set(pprob, 'cont', 'R', 0.1, 'PtMX', 300, 'indcs', [1 2 3]);
coco(pprob, 'torus1', [], 2, {'x' 'y' 'z'})

pprob = coco_add_func(pprob, 'angles', @angles, [], ...
  'zero', 'uidx', [1 2 3], 'u0', [-1; 1]);
pprob = coco_add_pars(pprob, 'phipar', [4 5], {'phi' 'theta'});

pprob = coco_set(pprob, 'cont', 'R', 0.1, 'PtMX', 300, 'indcs', [4 5]);
coco(pprob, 'torus2', [], 2, {'x' 'y' 'z' 'phi' 'theta'})

pprob = coco_set(pprob, 'cont', 'R', 0.1, 'PtMX', 300, 'indcs', [1 2 3]);
coco(pprob, 'torus3', [], 2, {'x' 'y' 'z' 'phi' 'theta'})

rmpath ../atlas2d_v2 % Remove atlas2d_v2 atlas algorithm from search path

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
coco_plot_bd('torus1', 'x', 'y', 'z')
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

% Projection onto angular variables
figure(2); clf; hold on
surf(X,Y,Z, 'EdgeColor', 0.7*[1 1 1], 'FaceColor', 0.8*[1 1 1]);
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'torus2', 'x', 'y', 'z')
coco_plot_bd('torus2', 'x', 'y', 'z')
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off

% Projection onto (x,y,z)-space with angles
figure(3); clf; hold on
surf(X,Y,Z, 'EdgeColor', 0.7*[1 1 1], 'FaceColor', 0.8*[1 1 1]);
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'torus3', 'x', 'y', 'z')
coco_plot_bd('torus3', 'x', 'y', 'z')
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
alpha 0.4; hold off
