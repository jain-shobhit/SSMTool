function figure_3_2
% Figure 3.2: A family of period-4 orbits of the Henon map obtained in
% Example 3.13. The family passes through a period-halving point marked
% with an open circle, where the period-4 family intersects a period-2
% family. A typical bifurcation diagram is shown in (a). Here, we graph a
% projection of the solution manifold onto the hyperplane with coordinates
% given by the parameter a and the x component of the first point on the
% orbit. The full state-space orbits are shown in (b).

% Note that the pairs (x1,y1) and (x3,y3), as well as the pairs (x2,y2) and
% (x4,y4), switch roles at the period-halving point.

% Generate data
if ~coco_exist('fig3.2', 'run')
  u0 = [1; 0.3; 1.275; -0.031; -0.656; 0.382; 0.952; ...
    -0.197; -0.103; 0.286; 1.275; -0.031];
  prob = period_B(u0, 4);
  prob = coco_add_pars(prob, 'pars', [1 2], {'a' 'b'});
  prob = coco_add_pars(prob, 'xy', 3:10, ...
    {'x1' 'y1' 'x2' 'y2' 'x3' 'y3' 'x4' 'y4'}, 'active');
  prob = coco_set(prob, 'cont', 'FP', true);
  coco(prob, 'fig3.2', [], 1, 'a', [0.8 1.2]);
end

% Extract data
bd  = coco_bd_read('fig3.2');
a   = coco_bd_col(bd, 'a');
x   = coco_bd_col(bd, 'x1');
lab = coco_bd_labs(bd, 'FP');
afp = coco_bd_val(bd, lab, 'a');
xfp = coco_bd_val(bd, lab, 'x1');

% Plot data: panel (a)
figure(1)
clf
hold on
box on
grid on
axis([0.9 1.21 0.7 1.3])

plot(a, x, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 15)
plot(afp, xfp, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', 'o', 'MarkerSize', 9.5, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

% Extract data
bd   = coco_bd_read('fig3.2');
x1   = coco_bd_col(bd, 'x1');
x2   = coco_bd_col(bd, 'x2');
y1   = coco_bd_col(bd, 'y1');
y2   = coco_bd_col(bd, 'y2');
lab  = coco_bd_labs(bd, 'FP');
x1fp = coco_bd_val(bd, lab, 'x1');
y1fp = coco_bd_val(bd, lab, 'y1');
x2fp = coco_bd_val(bd, lab, 'x2');
y2fp = coco_bd_val(bd, lab, 'y2');

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([-0.75 1.3 -0.25 0.4])

plot(x1, y1, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 15)
plot(x2, y2, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 15)
plot(x1fp, y1fp, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', 'o', 'MarkerSize', 9.5, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(x2fp, y2fp, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', 'o', 'MarkerSize', 9.5, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

end
