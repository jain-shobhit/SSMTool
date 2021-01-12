function figure_16_1
% Figure 16.1: Computation of the lemniscate curve in Example 16.1 with
% detection of atlas events. Two fold points at labels 6 and 10 and a
% branch point at label 3 are located. The covering closes at the end
% points 9 and 13.

% Generate data
if ~coco_exist('lemniscate', 'run')
  run demo_lemniscate
end

% Extract data
bd1 = coco_bd_read('lemniscate');
u  = coco_bd_col(bd1, {'x' 'y'});

% Plot data
figure(1)
clf
hold on
box on
grid on
axis([-1.25 1.25 -0.4 0.4])

plot(u(1,:), u(2,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd1, 'EP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd1, 'FP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd1, 'BP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

end
