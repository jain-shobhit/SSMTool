function figure_20_12
% Figure 20.12: Comparison of the high-period orbit (black) obtained in
% Fig. 20.1 with the low-period orbit (gray circles) used for constructing
% the initial solution guess in two different projectsion (a) and (b).
% Again, the number of mesh points along the excursion from the equilibrium
% is nealy identical for both orbits, although the mesh points move
% somewhat. The estimated error is 1.2506*10^(-5) for the low-period orbit
% and 1.5747*10^(-5) for the high-period orbit. Fig. 20.13 shows the
% adapted mesh and time profile of the high-period orbit.

% Note that machine-dependent round-off errors may affect the mesh
% discretization.

% Generate data
if ~(coco_exist('run1b', 'run') && coco_exist('run2b', 'run'))
  run demo_po_v3
end

coco_use_recipes_toolbox po_v3 coll_v6

% Extract data
sol0 = po_read_solution('', 'run1b', 6);
sol1 = po_read_solution('', 'run2b', 5);

[~, data0] = coll_read_solution('po.seg', 'run1b', 6);
[~, data1] = coll_read_solution('po.seg', 'run2b', 5);

idx1 = 1:4:data0.maps.NTST*4+1;
idx2 = 1:4:data1.maps.NTST*4+1;

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-0.3 0.25 -0.5 0.1])

plot(sol0.x(idx1,1), sol0.x(idx1,2), 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', [0.4 0.4 0.4], 'Marker', 'o', ...
  'MarkerSize', 8, 'MarkerFaceColor', 'white')
plot(sol1.x(:,1), sol1.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol1.x(idx2,1), sol1.x(idx2,2), 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', 'x', 'MarkerSize', 9.5) 

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([-0.5 0.1 -0.1 0.8])

plot(sol0.x(idx1,2), sol0.x(idx1,3), 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', [0.4 0.4 0.4], 'Marker', 'o', ...
  'MarkerSize', 8, 'MarkerFaceColor', 'white')
plot(sol1.x(:,2), sol1.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol1.x(idx2,2), sol1.x(idx2,3), 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', 'x', 'MarkerSize', 9.5) 

hold off

coco_use_recipes_toolbox

end
