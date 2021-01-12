% Examples of two-dimensional manifold atlases generated with the atlas_kd
% atlas algorithm.

%% Pillow
coco_func_data.pointers('set', [])
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);

% no computational boundaries
prob = coco_add_func(prob, 'pillow', @pillow, [], 'zero', ...
 'u0', [1; 0; 0]);
prob = coco_add_pars(prob, 'cartpar', 1:3, {'x', 'y', 'z'}, 'active');
prob = coco_set(prob, 'cont', 'PtMX', 5000);
coco(prob, 'pillow1', [], 2, {'x' 'y', 'z'});

% with computational boundaries
prob = coco_set(prob, 'cont', 'PtMX', 200, 'almax', 30);
coco(prob, 'pillow2', [], 2, {'x' 'y'}, {[0.5 1.5] [-0.5 0.5]});

atlas = coco_bd_read('pillow1', 'atlas');
figure(1); clf
plot_atlas_kd(atlas.charts, 3);
axis tight; axis equal;
view(3)
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [1.5 1.5 1.5], 'Style', 'local')

atlas = coco_bd_read('pillow2', 'atlas');
figure(2); clf
plot_atlas_kd(atlas.charts, 3);
axis tight; axis equal;
axis([0.5 1.5 -0.5 0.5 -inf inf])
view(3)
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [1.5 1.5 1.5], 'Style', 'local')

%% Ellipsoid
coco_func_data.pointers('set', [])
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);

% starting in the interior
prob = coco_add_func(prob, 'ellipsoid', @ellipsoid, [], 'zero', ...
  'u0', [1.3; 0.7; 0.4]);
prob = coco_add_pars(prob, 'cartpar', 1:3, {'x' 'y' 'z'}, 'active');
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'almax', 20);
coco(prob, 'ellipsoid1', [], 2, {'x' 'y' 'z'}, {[] [0 1] [0 1]});

% starting on the boundary
prob = coco_add_func(prob, 'ellipsoid', @ellipsoid, [], 'zero', ...
  'x0', [1.49749; 0.1; 0]);
prob = coco_add_pars(prob, 'cartpar', 1:3, {'x' 'y' 'z'}, 'active');
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'almax', 10);
coco(prob, 'ellipsoid2', [], 2, {'x' 'y' 'z'}, {[] [0 1] [0 1]});

atlas = coco_bd_read('ellipsoid1', 'atlas');
figure(3); clf
plot_atlas_kd(atlas.charts, 3);
axis tight; axis equal; view(-60,20)
axis([-inf inf 0 1 0 1])
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [1.5 1.5 1.5], 'Style', 'local')

atlas = coco_bd_read('ellipsoid2', 'atlas');
figure(4); clf
plot_atlas_kd(atlas.charts, 3);
axis tight; axis equal; view(-60,20)
axis([-inf inf 0 1 0 1])
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [1.5 1.5 1.5], 'Style', 'local')

%% Cylinder
coco_func_data.pointers('set', [])
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);

prob = coco_add_func(prob, 'cylinder', @cylinder, [], ...
  'zero', 'u0', [1; 0; 0] + sqrt([0.5; 0.55; 0]));
prob = coco_add_func(prob, 'wedge', @wedge, [], 'active', 'phi', ...
  'uidx', 1:3);
prob = coco_add_pars(prob, 'cartpar', 1:3, {'x' 'y' 'z'}, 'active');
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'almax', 20);
coco(prob, 'cylinder', [], 2, 'phi', [-0.5 0.5]);

atlas = coco_bd_read('cylinder', 'atlas');
figure(5); clf
plot_atlas_kd(atlas.charts, 2, 3, 4);
axis tight; axis equal; view(30,10)
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [2.5 0 0], 'Style', 'local')

%% Sphere
coco_func_data.pointers('set', [])
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);

prob = coco_add_func(prob, 'sphere', @sphere, [], 'zero', ...
  'u0', [1; 1; 1]);
prob = coco_add_pars(prob, 'cartpar', 1:3, {'x', 'y', 'z'}, 'active');
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'almax', 20);
prob = coco_add_event(prob, 'UZ', 'z', 0.5);
coco(prob, 'sphere', [], 2, 'x', [0.5 1.5]);

[atlas, bd] = coco_bd_read('sphere', 'atlas', 'bd');
figure(6); clf;
plot_atlas_kd(atlas.charts, 3);
hold on
idx = coco_bd_idxs(bd, 'UZ');
x   = coco_bd_col(bd, {'x' 'y' 'z'});
plot3(x(1,idx), x(2,idx), x(3,idx), 'b.', 'MarkerSize', 21);
hold off; axis equal; axis tight; view(60,30)
axis([0.5 1.5 -inf inf -inf inf])
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [1.5 0 0], 'Style', 'local')

%% Klein bottle
coco_func_data.pointers('set', [])
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);

prob = coco_add_func(prob, 'klein', @klein, [], 'zero', ...
  'u0', [0; -1-sqrt(2); 0]);
prob = coco_add_pars(prob, 'cartpar', 1:3, {'x' 'y' 'z'}, 'active');
prob = coco_set(prob, 'cont', 'PtMX', 2000, 'almax', 40, 'MaxRes', 20);
coco(prob, 'klein', [], 2);

atlas = coco_bd_read('klein', 'atlas');
figure(7); clf
plot_atlas_kd(atlas.charts, 3);
axis equal; axis tight; axis off; view(-50,4)
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [1.5 0 0], 'Style', 'local')

%% Torus
coco_func_data.pointers('set', [])
prob = coco_prob();
prob = coco_set(prob, 'all', 'CleanData', true);

prob = coco_add_func(prob, 'torus', @torus, [], 'zero', 'u0', [0; 2; 1]);
prob = coco_set(prob, 'cont', 'PtMX', 1000, 'almax', 40, 'R_max', 0.4);

prob = coco_add_pars(prob, 'cartpar', [1 2 3], {'x' 'y' 'z'});
coco(prob, 'torus1', [], 2, {'x' 'y' 'z'})

prob = coco_add_func(prob, 'angles', @angles, [], ...
  'zero', 'uidx', [1 2 3], 'u0', [pi/2; pi/2]);
prob = coco_add_pars(prob, 'phipar', [4 5], {'theta' 'phi'});
coco(prob, 'torus2', [], 2, {'theta' 'phi'})

prob = coco_add_func(prob, 'angles', @angles, [], ...
  'zero', 'uidx', [1 2 3], 'u0', [pi/2; pi/2]);
prob = coco_add_pars(prob, 'phipar', [3 4 5], {'z' 'theta' 'phi'});
coco(prob, 'torus3', [], 2, {'z' 'theta' 'phi'}, [0 2])

n = 20;
angle = 0:2*pi/n:2*pi;
X = (2+0.99.*cos(angle)')*cos(angle);
Y = (2+0.99.*cos(angle)')*sin(angle);
Z = 0.99.*sin(angle)'*ones(1,n+1);

% Projection onto (x,y,z)-space without angles.
figure(8); clf; hold on
surf(X,Y,Z, 'EdgeColor', 0.7*[1 1 1], 'FaceColor', 0.8*[1 1 1]);
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'torus1', 'x', 'y', 'z')
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [3 2 2], 'Style', 'local')

figure(9); clf; hold on
surf(X,Y,Z, 'EdgeColor', 0.7*[1 1 1], 'FaceColor', 0.8*[1 1 1]);
atlas = coco_bd_read('torus1', 'atlas');
plot_atlas_kd(atlas.charts, 3)
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [3 2 2], 'Style', 'local')

% Projection onto angular variables
figure(10); clf; hold on
surf(X,Y,Z, 'EdgeColor', 0.7*[1 1 1], 'FaceColor', 0.8*[1 1 1]);
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'torus3', ...
  {'theta' 'phi'}, @(phi,theta) (2+cos(theta)).*cos(phi), ...
  {'theta' 'phi'}, @(phi,theta) (2+cos(theta)).*sin(phi), 'z')
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [3 2 2], 'Style', 'local')

figure(11); clf; hold on
atlas = coco_bd_read('torus2', 'atlas');
plot_atlas_kd(atlas.charts, 2)
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
grid on
box on
set(gca,'LineWidth',2)
hold off

% Projection onto (x,y,z)-space with angles
figure(12); clf; hold on
surf(X,Y,Z, 'EdgeColor', 0.7*[1 1 1], 'FaceColor', 0.8*[1 1 1]);
thm = struc([]);
thm.lspec = {'ro', 'LineWidth', 1};
coco_plot_bd(thm, 'torus3', ...
  {'theta' 'phi'}, @(th,phi) (2+cos(phi)).*cos(th), ...
  {'theta' 'phi'}, @(th,phi) (2+cos(phi)).*sin(th), ...
  'z')
axis([-3 3 -3 3 -2 2]); axis equal; grid on; view(60,30); 
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [3 2 2], 'Style', 'local')

atlas = coco_bd_read('torus3', 'atlas');
figure(13); clf; hold on
plot_atlas_kd(atlas.charts, 3)
axis([-3 3 -3 3 0 2]); axis equal; grid on; view(60,30); 
grid on
box on
set(gca,'LineWidth',2)
hold off
alpha 0.4
light('Position', [3 2 2], 'Style', 'local')
