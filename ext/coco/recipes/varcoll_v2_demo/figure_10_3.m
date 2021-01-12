function figure_10_3
% Figure 10.3: Construction of an initial approximation of an orbit
% connecting an equilibrium and a periodic orbit in the Lorenz system given
% by the vector field in Eq. (10.87) following the homotopy approach
% described in Sect. 10.2.2. A state-space representation of the
% three-segment solution, consisting of a periodic orbit (gray) and two
% zero-length segments (black dots), that is used to initialize Stage II of
% the homotopy is shown in panel (a), together with the hyperplane Sigma
% that separates the periodic orbit from the equilibrium. In Stage II, we
% grow an orbit in W_0^u until it terminates on Sigma, as shown in (b). In
% the subsequent Stage III we grow an orbit in W_per^s in a similar way;
% see Fig. 10.4.

% Generate data
if ~(coco_exist('runHopf', 'run') && coco_exist('run1', 'run') ...
    && coco_exist('run2', 'run') && coco_exist('run3', 'run') ...
    && coco_exist('run4', 'run') && coco_exist('run5', 'run')) 
  run demo_lorenz
end

coco_use_recipes_toolbox coll_v2

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
axis([-10 40 -20 30 -10 40])
view([35 25])
camproj('perspective')

sol = coll_read_solution('po.seg', 'run1', 1); % Extract periodic orbit segment
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)
sol = coll_read_solution('col1', 'run1', 1); % Extract segment in unstable manifold of equilibrium
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
sol = coll_read_solution('col2', 'run1', 1); % Extract segment in stable manifold of orbit
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

[X Y] = meshgrid([-10 40], [-20 30]);
Z     = @(n,x0, x,y) (x0'*n)/n(3) - (n(1)/n(3))*x - (n(2)/n(3))*y;
n     = [0;-1;1];
x0    = [20;20;30];
surf(X,Y,Z(n,x0,X,Y), 'facecolor', 0.5*[1 1 1], 'facealpha', ...
  0.5, 'linestyle', '-')

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
axis([-10 40 -20 30 -10 40])
view([35 25])
camproj('perspective')

sol = coll_read_solution('po.seg', 'run1', 6); % Extract periodic orbit segment
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)
sol = coll_read_solution('col1', 'run1', 6); % Extract segment in unstable manifold of equilibrium
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
sol = coll_read_solution('col2', 'run1', 6); % Extract segment in stable manifold of orbit
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

[X Y] = meshgrid([-10 40], [-20 30]);
Z     = @(n,x0, x,y) (x0'*n)/n(3) - (n(1)/n(3))*x - (n(2)/n(3))*y;
n     = [0;-1;1];
x0    = [20;20;30];
surf(X,Y,Z(n,x0,X,Y), 'facecolor', 0.5*[1 1 1], 'facealpha', ...
  0.5, 'linestyle', '-')

hold off

coco_use_recipes_toolbox

end
