function figure_10_4
% Figure 10.4: Panel (a) shows the three-segment solution after completing
% Stage III, i.e., growing an orbit in W_per^s starting at the solution
% shown in Fig. 10.3(b). Here, the end point of the orbit segment in W_0^u
% and the starting point of the orbit segment in W_per^s both terminate on
% Sigma. Although W_per^s is 2-dimensional, the connecting orbit is unique.
% To obtain an initial approximation of the connecting orbit, we first
% sweep W_per^s in Stage IV and compute a set of orbit segments that cover
% the manifold sufficiently densely (b). From this family of orbits we
% select the one that terminates closest to the end point of the segment in
% W_0^u; see Fig. 10.5.

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
axis([-10 40 -20 30 -10 50])
view([35 25])
camproj('perspective')

sol = coll_read_solution('po.seg', 'run2', 9); % Extract periodic orbit segment
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)
sol = coll_read_solution('col1', 'run2', 9); % Extract segment in unstable manifold of equilibrium
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
sol = coll_read_solution('col2', 'run2', 9); % Extract segment in stable manifold of orbit
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
axis([-10 40 -20 30 -10 50])
view([35 25])
camproj('perspective')

XX = [];
YY = [];
ZZ = [];

bd = coco_bd_read('run3');
for lab = coco_bd_labs(bd)
  sol = coll_read_solution('col2', 'run3', lab);
  XX = [XX sol.x(:,1)]; %#ok<AGROW>
  YY = [YY sol.x(:,2)]; %#ok<AGROW>
  ZZ = [ZZ sol.x(:,3)]; %#ok<AGROW>
end

sol = coll_read_solution('po.seg', 'run3', 1); % Extract periodic orbit segment
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)
sol = coll_read_solution('col1', 'run3', 1); % Extract segment in unstable manifold of equilibrium
plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

surf(XX,YY,ZZ, 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 1, ...
  'MeshStyle', 'column', 'LineStyle', '-', 'EdgeColor', 0.6*[1 1 1], ...
  'LineWidth', 0.5)

plot3(XX(1,:), YY(1,:), ZZ(1,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'white')

[X Y] = meshgrid([-10 40], [-20 30]);
Z = @(n,x0, x,y) (x0'*n)/n(3) - (n(1)/n(3))*x - (n(2)/n(3))*y;
n = [0;-3;3];
x0 = [20;20;30];
surf(X,Y,Z(n,x0,X,Y), 'facecolor', 0.5*[1 1 1], 'facealpha', 0.5, ...
  'linestyle', '-')

hold off

coco_use_recipes_toolbox

end
