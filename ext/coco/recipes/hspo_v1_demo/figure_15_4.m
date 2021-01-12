function figure_15_4
% Figure 15.4: Bifurcation diagram of the impact oscillator from Example
% 15.2 with ||U||:=||(u,mu_J)||. Continuation starts with an impacting
% orbit at label 1 and passes through a grazing point. The part of the
% curve between labels 6 and 11 corresponds to nonphysical orbits
% penetrating the impact surface. Some representative orbits are shown in
% Fig. 15.5.

% Generate data
if ~(coco_exist('impact1', 'run') && coco_exist('impact2', 'run'))
  run demo_impact
end

% Extract data
bd1 = coco_bd_read('impact1');
A   = coco_bd_col(bd1, 'A');
y   = coco_bd_col(bd1, '||U||');

% Plot data
figure(1)
clf
hold on
grid on
box on
axis([0.05 1.05 21 31])

plot(A, y, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd1, 'EP');
plot(A(idx), y(idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd1, 'GR');
plot(A(idx), y(idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

end
