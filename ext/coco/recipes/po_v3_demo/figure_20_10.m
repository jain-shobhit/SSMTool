function figure_20_10
% Figure 20.10: Time profiles of the high-period orbit constructed in Fig.
% 20.9 on a moving mesh shown over the full period (a) and in two
% subsequent zooms into the excursion from the equilibrium (b) and (c). The
% distribution of mesh points is indicated with circles at the bottom.

% Note that machine-dependent round-off errors may affect the mesh
% discretization.

% Generate data
if ~(coco_exist('run1a', 'run') && coco_exist('run2a', 'run'))
  run demo_po_v3
end

coco_use_recipes_toolbox po_v3 coll_v5

% Extract data
sol1 = po_read_solution('', 'run2a', 2);

N   = numel(sol1.t);
N2  = 4*round((N-1)/8)+1;
idx = [N2:N 2:N2];
tt  = [sol1.t(N2:N)-sol1.t(N2) ; sol1.t(N)+sol1.t(2:N2)-sol1.t(N2)];
tt  = tt/tt(end);
idx2 = 1:4:120*4+1;

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0 1 -0.7 0.8])

plot(tt, sol1.x(idx,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(tt(idx2), 0*sol1.x(idx2,1)-0.6, 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0.474 0.485 -0.7 0.8])

plot(tt, sol1.x(idx,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(tt(idx2), 0*sol1.x(idx2,1)-0.6, 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (c)
figure(3)
clf
hold on
grid on
box on
axis([0.4801 0.4817 -0.7 0.8])

plot(tt, sol1.x(idx,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(tt(idx2), 0*sol1.x(idx2,1)-0.6, 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

coco_use_recipes_toolbox

end
