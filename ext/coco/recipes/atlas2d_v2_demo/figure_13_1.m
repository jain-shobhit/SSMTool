function figure_13_1
% Figure 13.1: Coverings of a unit sphere obtained in Example 13.1 with the
% modified atlas_2d_min atlas algorithm from Sect. 13.2. The first run
% results in an incomplete covering shown in panel (a): the algorithm
% terminates prematurely after covering a strip around the equator. The
% second run, with a different choice for the number of continuation
% directions for each chart, results in a complete cover (b). The dots mark
% the position of the base points of the charts in the final atlas.

% Generate data
if ~coco_exist('sphere1', 'run') || ~coco_exist('sphere2', 'run')
  run demo_atlas2d_v2
end

d = cd('..');
N = 16;
M = 5;
[X Y Z]=sphere(N*M);
idx = [1 (1:N)*M+1];
cd(d);

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
axis([0 2 -1 1 -1 1], 'equal')
view([60 30])

bd = coco_bd_read('sphere1');        % Extract bifurcation data
x  = coco_bd_col(bd, {'x' 'y' 'z'}); % Extract column data
plot3(x(1,:), x(2,:), x(3,:), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black')
surf(0.981*X(idx,:)+1,0.981*Y(idx,:),0.981*Z(idx,:), 'FaceColor', 'none', ...
  'EdgeColor', 0.5*[1 1 1], 'MeshStyle', 'row', 'LineWidth', 1)
surf(0.981*X(:,idx)+1,0.981*Y(:,idx),0.981*Z(:,idx), 'FaceColor', 'none', ...
  'EdgeColor', 0.5*[1 1 1], 'MeshStyle', 'column', 'LineWidth', 1)
surf(0.98*X+1,0.98*Y,0.98*Z, 'EdgeColor', 'none', ...
  'FaceColor', 0.8*[1 1 1], 'FaceAlpha', 0.999)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
axis([0 2 -1 1 -1 1], 'equal')
view([60 30])

bd = coco_bd_read('sphere2');        % Extract bifurcation data
x  = coco_bd_col(bd, {'x' 'y' 'z'}); % Extract column data
plot3(x(1,:), x(2,:), x(3,:), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black')
surf(0.981*X(idx,:)+1,0.981*Y(idx,:),0.981*Z(idx,:), 'FaceColor', 'none', ...
  'EdgeColor', 0.5*[1 1 1], 'MeshStyle', 'row', 'LineWidth', 1)
surf(0.981*X(:,idx)+1,0.981*Y(:,idx),0.981*Z(:,idx), 'FaceColor', 'none', ...
  'EdgeColor', 0.5*[1 1 1], 'MeshStyle', 'column', 'LineWidth', 1)
surf(0.98*X+1,0.98*Y,0.98*Z, 'EdgeColor', 'none', ...
  'FaceColor', 0.8*[1 1 1], 'FaceAlpha', 0.999)

hold off

end
