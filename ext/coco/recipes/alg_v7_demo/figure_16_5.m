function figure_16_5
% Figure 16.5: The covering of the cusp surface obtained in the last run of
% Sect. 16.2.3; fold points are marked as black dots. The 2-dimensional
% atlas algorithm locates a fold point whenever a curve segment intersects
% an event surface transversally. In panel (a), one can see that events
% between neighboring charts may remain undetected if no connecting curve
% segment was constructed by the atlas algorithm. Projecting the manifold
% onto the (kappa,lambda) parameter plane again reveals the fold curve as
% the border curve between different shades of gray (b).

% Generate data
if ~(coco_exist('cuspsurface', 'run'))
  run demo_alg_v7
end
coco_use_recipes_toolbox atlas2d_v6

% Extract data
[atlas bd] = coco_bd_read('cuspsurface', 'atlas', 'bd');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
axis([-1 1 -0.5 0.5 0 1.2])
view([100 25])

[tri X] = plot_trisurf(atlas.charts, 3, 2, 1);
x0      = repmat([0 0 0], size(X,1), 1);
n       = [1 1 1];
C       = n*(X-x0)';
cmap    = repmat(linspace(0.55,0.95,100)', 1, 3);
colormap(cmap);
surf([0 0;0 0], [0 0;0 0], [0 0;0 0], 'FaceColor', 'white', ...
  'EdgeColor', 'white', 'FaceAlpha', 0)
trisurf(tri, X(:,1), X(:,2), abs(X(:,3)), C, 'FaceColor', 'interp', ...
  'EdgeColor', 0.5*[1 1 1])

x  = coco_bd_col(bd, '||alg.x||');
ka = coco_bd_col(bd, 'ka');
la = coco_bd_col(bd, 'la');
idx = coco_bd_idxs(bd, 'FO');
plot3(la(idx)+0.05, ka(idx), x(idx), 'ko', 'MarkerSize', 6, ...
  'MarkerFaceColor', 'k')

axis tight
hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
view([0 90])

surf([0 0;0 0], [0 0;0 0], [0 0;0 0], 'FaceColor', 'white', ...
  'EdgeColor', 'white', 'FaceAlpha', 0)
trisurf(tri, X(:,2), X(:,1), X(:,3), 'FaceColor', 0.8*[1 1 1], ...
  'EdgeColor', 0.7*[1 1 1], 'LineWidth', 0.5, 'FaceAlpha', 0.35)
x  = coco_bd_col(bd, '||alg.x||');
ka = coco_bd_col(bd, 'ka');
la = coco_bd_col(bd, 'la');
idx = coco_bd_idxs(bd, 'FO');
plot3(ka(idx), la(idx), x(idx)+1.3, 'ko', 'MarkerSize', 6, ...
  'MarkerFaceColor', 'k')

axis tight
hold off

coco_use_recipes_toolbox

end
