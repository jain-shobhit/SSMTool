function figure_20_13
% Figure 20.13: Time profiles of the high-period orbit constructed in Fig.
% 20.11 on a moving mesh with variable discretization order shown over the
% full period (a) and in two subsequent zooms into the excursion from the
% equilibrium (b) and (c). The distribution of mesh points is indicated
% with circles at the bottom. A comparison between panel (a) and Fig.
% 20.10(a) indicates that the difference in discretization order of 41 mesh
% intervals is due mainly to a coarser mesh along the part with
% near-constant dynamics. The ratio between the largest and smallest mesh
% intervals, max(kappa)/min(kappa)=9245, is here almost twice as large as
% for the moving mesh with fixed discretization order used in Fig. 20.10.

% Note that machine-dependent round-off errors may affect the mesh
% discretization.

% Generate data
if ~(coco_exist('run1b', 'run') && coco_exist('run2b', 'run'))
  run demo_po_v3
end

coco_use_recipes_toolbox po_v3 coll_v6

% Extract data
sol1 = po_read_solution('', 'run2b', 5);

[~, data] = coll_read_solution('po.seg', 'run2b', 5);
% fprintf('max(data.mesh.ka)/min(data.mesh.ka) = %.3e\n', ...
%    max(data.mesh.ka)/min(data.mesh.ka));

N   = numel(sol1.t);
N2  = 4*round((N-1)/8)+1;
idx = [N2:N 2:N2];
tt  = [sol1.t(N2:N)-sol1.t(N2) ; sol1.t(N)+sol1.t(2:N2)-sol1.t(N2)];
tt  = tt/tt(end);
idx2 = 1:4:data.maps.NTST*4+1;

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
axis([0.602 0.613 -0.7 0.8])

plot(tt, sol1.x(idx,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(tt(idx2), 0*sol1.x(idx2,1)-0.6, 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (c)
figure(3)
clf
hold on
grid on
box on
axis([0.6075 0.6091 -0.7 0.8])

plot(tt, sol1.x(idx,:), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(tt(idx2), 0*sol1.x(idx2,1)-0.6, 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

coco_use_recipes_toolbox

end
