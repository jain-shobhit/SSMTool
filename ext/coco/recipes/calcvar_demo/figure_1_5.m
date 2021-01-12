function figure_1_5
% Figure 1.5: A branch of solutions of Eqs. (1.12-1.13) obtained using
% numerical continuation as described in Sect. 1.3.1. Here, dots indicate
% computed points, which are connected with straight-line segments. The
% circle denotes the initial point at a=1, b=0, and Y=cosh 1. Compare the
% numerical results with the analysis shown in Fig. 1.2.

% Generate data
if ~coco_exist('run_alg', 'run')
  run demo_algebraic
end

% Extract data
bd = coco_bd_read('run_alg');
Y  = coco_bd_col(bd, 'Y');
a  = coco_bd_col(bd, 'a');
b  = coco_bd_col(bd, 'b');

Y0 = coco_bd_val(bd, 1, 'Y');
a0 = coco_bd_val(bd, 1, 'a');
b0 = coco_bd_val(bd, 1, 'b');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0 10 0 8])

plot(Y, a, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 15)
plot(Y0, a0, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 10 -0.7 0.7])

plot(Y, b, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 15)
plot(Y0, b0, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

end
