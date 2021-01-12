function figure_20_6
% Figure 20.6: Continuation of periodic orbits of the dynamical system
% given by the vector field in Eq. (18.5) using an equidistributed moving
% mesh with adaptive changes to the discretization order, as described in
% Sect. 20.2.2. Here, a single continuation run suffices to compute entire
% solution family. Due to the conservative definition of the adaptation
% window, the algorithm uses about 50% more mesh points than the analysis
% shown in Fig. 20.5. Panel (a) shows sample orbits obtained during
% continuation, and panels (b) and (c) illustrate the effects of
% adaptation.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run'))
  run demo_coll_v6
end

coco_use_recipes_toolbox po_v3 coll_v6

% Extract data
bd1 = coco_bd_read('run1'); % Extract bifurcation data
bd2 = coco_bd_read('run2'); % Extract bifurcation data

sol1 = po_read_solution('', 'run1', 3);
[~, data1] = coll_read_solution('po.seg', 'run1', 3);
sol2 = po_read_solution('', 'run1', 7);
[~, data2] = coll_read_solution('po.seg', 'run1', 7);
sol3 = po_read_solution('', 'run1', 10);
[~, data3] = coll_read_solution('po.seg', 'run1', 10);

coll = data1.coll;
TOL1 = coll.TOLINC;
TOL2 = coll.TOLDEC;

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-3.0833 3.0833 -1.0794*(1-1/6) 1.0794*(1-1/6)]*1.05)

idx = 1:4:data3.maps.NTST*4+1;
plot(sol3.x(idx,1), sol3.x(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(sol3.x(:,1), sol3.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

idx = 1:4:data2.maps.NTST*4+1;
plot(sol2.x(idx,1), sol2.x(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(sol2.x(:,1), sol2.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

idx = 1:4:data1.maps.NTST*4+1;
plot(sol1.x(idx,1), sol1.x(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(sol1.x(:,1), sol1.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 20 0 coll.TOL])

p1 = coco_bd_col(bd1, 'eps');
err1 = coco_bd_col(bd1, 'po.seg.coll.err'); % Extract column data
N1 = coco_bd_col(bd1, 'po.seg.coll.NTST');  % Extract column data

p2 = coco_bd_col(bd2, 'eps'); % Extract column data
err2 = coco_bd_col(bd2, 'po.seg.coll.err'); % Extract column data
N2 = coco_bd_col(bd2, 'po.seg.coll.NTST');  % Extract column data

plot([0 20], [TOL1 TOL1], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', [0.3 0.3 0.3])
plot([0 20], [TOL2 TOL2], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', [0.3 0.3 0.3])

plot(p1, err1, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)
plot(p2, err2, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (c)
figure(3)
clf
hold on
grid on
box on
axis([0 20 0 40])

plot(p1, N1, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)
plot(p2, N2, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)

hold off

coco_use_recipes_toolbox

end
