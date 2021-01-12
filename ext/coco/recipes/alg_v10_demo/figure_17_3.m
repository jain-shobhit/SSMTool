function figure_17_3
% Figure 17.3: Repeating the computation of Example 17.1 using an event
% handler, as shown in Sect. 17.1.3, allows for a distinction between Hopf
% and neutral saddle point; compare with Fig. 17.1(b). Panel (a) shows only
% Hopf and fold points, while panel (b) also shows neutral saddle points
% marked with x. The gray curve is the locus of zeros of psi_HB given by
% Eq. (17.4).

% Generate data
if ~(coco_exist('run', 'run'))
  run demo_alg_v10
end

% Extract data
bd = coco_bd_read('run', 'bd');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0 0.5 0 0.255])

p1 = coco_bd_col(bd, 'p1');
p2 = coco_bd_col(bd, 'p2');
idx1 = coco_bd_idxs(bd, 'HB');
idx2 = coco_bd_idxs(bd, 'FO');
idx3 = coco_bd_idxs(bd, 'NSad');

x = linspace(0,0.25,100);
y = (1+sqrt(1-4*x)+2*x)./(4+2*x);
plot(y, x, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
y = (1-sqrt(1-4*x)+2*x)./(4+2*x);
plot(y, x, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])

plot(p1(idx1), p2(idx1), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'black')
plot(p1(idx2), p2(idx2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 9.5, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 0.5 0 0.255])

x = linspace(0,0.25,100);
y = (1+sqrt(1-4*x)+2*x)./(4+2*x);
plot(y, x, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
y = (1-sqrt(1-4*x)+2*x)./(4+2*x);
plot(y, x, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])

plot(p1(idx1),p2(idx1), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'black')
plot(p1(idx3),p2(idx3), 'LineStyle', 'none', 'LineWidth', sqrt(2), ...
  'Color', 'black', 'Marker', 'x', 'MarkerSize', 9.5)
plot(p1(idx2),p2(idx2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 9.5, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

end
