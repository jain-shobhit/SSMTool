function figure_20_8
% Figure 20.8: Time profiles of the high-period orbit constructed in Fig.
% 20.7 on a uniform mesh shown over the full period (a) and in a zoom into
% the excursion from the equilibrium (b). The distribution of mesh points
% is indicated with circles at the bottom. In panel (a), the mesh is so
% dense that it is not possible to see individual mesh intervals.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run'))
  run demo_po_v3
end

coco_use_recipes_toolbox po_v1 coll_v3

% Extract data
sol1 = po_read_solution('', 'run2', 1);

N   = numel(sol1.t);
N2  = round(N/2);
idx = [N2:N 2:N2];
tt  = [sol1.t(N2:N)-sol1.t(N2) ; sol1.t(N)+sol1.t(2:N2)-sol1.t(N2)];
tt  = tt/tt(end);
idx2 = 1:4:1250*4+1;

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
  'LineWidth', 2, 'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0.45 0.55 -0.7 0.8])

plot(tt, sol1.x(idx,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(tt(idx2), 0*sol1.x(idx2,1)-0.6, 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

coco_use_recipes_toolbox

end
