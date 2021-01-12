function figure_15_6
% Figure 15.6: Grazing curve in the (A,omega) parameter plane computed for
% the impact oscillator in Example 15.2. For the system under
% consideration, this curve subdivides the parameter plane locally into a
% region of parameter values for which nonimpacting orbits exist (left-hand
% side) from parameter values for which orbits have at least one impact
% (right-hand side). Some representive grazing orbits are shown in Fig.
% 15.7.

% Generate data
if ~(coco_exist('impact1', 'run') && coco_exist('impact2', 'run'))
  run demo_impact
end

% Extract data
bd1 = coco_bd_read('impact2');
A   = coco_bd_col(bd1, 'A');
y   = coco_bd_col(bd1, 'w');

% Plot data
figure(1)
clf
hold on
grid on
box on
axis([0 1 0 1.5])

plot(A, y, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd1, 'EP');
plot(A(idx), y(idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

end
