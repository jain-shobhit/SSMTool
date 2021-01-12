function figure_16_4
% Figure 16.4: Fold-point detection and continuation for the cusp normal
% form considered in Sect. 16.2.3 using the monitor function defined in
% Sect. 16.2.1 embedded in the 'alg' toolbox. Panels (a) and (b) show the
% result of computing the solution to the restricted continuation problem
% illustrated in Fig. 15.2 with the monitor function added as either
% regular or active, respectively, and with ||U||:=||(u,mu_J)||. The
% difference in appearance is due to the exclusion or inclusion of the
% value of the monitor function in the vector of continuation parameters
% mu. Two fold points are located in each case. Restarting at one of these
% fold points while restricting the monitor function to 0 results in a
% covering of the locus of fold points (c); see also Fig. 15.1.

% Generate data
if ~(coco_exist('regular', 'run') && coco_exist('active', 'run') ...
    && coco_exist('cusp', 'run'))
  run demo_alg_v7
end

% Extract data
bd1 = coco_bd_read('regular');
bd2 = coco_bd_read('active');
bd3 = coco_bd_read('cusp');

% Plot data: panel (a)
figure(1)
clf
hold on
box on
grid on
axis([-0.5 0.5 0.4 1.5])

x  = coco_bd_col(bd1, '||U||'); % Extract column data
ka = coco_bd_col(bd1, 'ka');    % Extract column data
plot(ka, x, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd1, 'EP'); % Extract row indices
plot(ka(idx), x(idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd1, 'FO'); % Extract row indices
plot(ka(idx), x(idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
box on
grid on
axis([-0.5 0.5 0.5 3])

x  = coco_bd_col(bd2, '||U||'); % Extract column data
ka = coco_bd_col(bd2, 'ka');    % Extract column data
plot(ka, x, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd2, 'EP'); % Extract row indices
plot(ka(idx), x(idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')
idx = coco_bd_idxs(bd2, 'FO'); % Extract row indices
plot(ka(idx), x(idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (c)
figure(3)
clf
hold on
box on
grid on
axis([-0.5 0.5 0 1.3])

la = coco_bd_col(bd3, 'la'); % Extract column data
ka = coco_bd_col(bd3, 'ka'); % Extract column data
plot(ka, la, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
idx = coco_bd_idxs(bd3, 'EP'); % Extract row indices
plot(ka(idx), la(idx), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

end
