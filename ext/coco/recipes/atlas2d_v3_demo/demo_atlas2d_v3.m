coco_use_recipes_toolbox atlas2d_v3 % Add atlas2d_v3 atlas algorithm to search path

%% Analogous to Example 13.1 [code not in Recipes for Continuation]

% The continuation problem structure encoded below corresponds to a single
% zero function in terms of three continuation variables, a family of three
% monitor functions that evaluate to the continuation variables, and the
% corresponding inactive continuation parameters 'x', 'y', and 'z'. Its
% dimensional deficit is -1. A two-dimensional solution manifold results by
% releasing 'x', 'y', and 'z' and allowing these to vary during
% continuation.

% In the first run, the number of available directions for
% continuation from individual charts is 6, but the algorithm terminates
% prematurely after going once around a great circle. In the second run,
% the number of available directions for continuation is only 4. The
% algorithm covers the entire sphere without redundancy.

prob = coco_add_func(coco_prob(), 'sphere', @sphere, [], ...
  'zero', 'u0', [2;0;0]);
prob = coco_add_pars(prob, '', [1 2 3], {'x' 'y' 'z'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'h', 0.5, 'almax', 35); % Number of available directions = 6
coco(prob, 'sphere1', [], 2, {'x' 'y' 'z'});

prob = coco_set(prob, 'cont', 'Ndirs', 4, 'PtMX', 200); % Number of available directions = 4
coco(prob, 'sphere2', [], 2, {'x' 'y' 'z'});

figure(1)
clf
p = cd('../');
[X Y Z] = sphere(16); % Execute Matlab function sphere, rather than local function
cd(p);
bd = coco_bd_read('sphere1');        % Extract bifurcation data
x  = coco_bd_col(bd, {'x' 'y' 'z'}); % Extract column data
subplot(1,2,1)
plot3(x(1,:), x(2,:), x(3,:), 'k.', 'MarkerSize', 12);
hold on
surf(0.99*X+1,0.99*Y,0.99*Z, 'EdgeColor', 0.7*[1 1 1], ...
  'FaceColor', 0.8*[1 1 1]);
hold off
view(60,30)
grid on
axis([0 2 -1 1 -1 1]);
axis equal
drawnow

bd = coco_bd_read('sphere2');        % Extract bifurcation data
x  = coco_bd_col(bd, {'x' 'y' 'z'}); % Extract column data
subplot(1,2,2)
plot3(x(1,:), x(2,:), x(3,:), 'k.', 'MarkerSize', 12);
hold on
surf(0.99*X+1,0.99*Y,0.99*Z, 'EdgeColor', 0.7*[1 1 1], ...
  'FaceColor', 0.8*[1 1 1]);
hold off
view(60,30)
grid on
axis([0 2 -1 1 -1 1]);
axis equal
drawnow

coco_use_recipes_toolbox % Remove atlas2d_v3 atlas algorithm from search path
