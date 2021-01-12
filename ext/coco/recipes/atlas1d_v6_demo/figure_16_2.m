function figure_16_2
% Figure 16.2: Frequency response curve of the Duffing oscillator given by
% the vector field in Eq. (16.10) for A=2.5, lambda=0.2, and
% alpha=epsilon=1, as computed in Example 16.2 with ||U||:=||(u,mu_J)||. A
% number of atlas events are detected along this curve, as shown in panel
% (a) and the enlargements (b) and (c). For sample solutions at labels 16,
% 35, and 58 along the curve we observe a characteristic phase shift, as
% illustrated in Fig. 16.3.

% Generate data
if ~coco_exist('duffing', 'run')
  run demo_duffing
end

% Extract data
bd2 = coco_bd_read('duffing');
u  = coco_bd_col(bd2, {'om' '||U||'});

% Plot data: panel (a)
figure(1)
clf
hold on
box on
grid on
axis([0.5 3.5 12 120])

plot(u(1,:), u(2,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd2, 'EP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd2, 'FP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd2, 'BP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
box on
grid on
axis([0.57 1.08 22.55 27.3])

plot(u(1,:), u(2,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd2, 'EP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd2, 'FP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd2, 'BP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (c)
figure(3)
clf
hold on
box on
grid on
axis([0.588 0.6 25.9 26.9])

plot(u(1,:), u(2,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd2, 'EP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd2, 'FP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd2, 'BP');
plot(u(1,idx),u(2,idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

end
