function figure_18_5
% Figure 18.5: Continuation of periodic orbits of the dynamical system
% given by the vector field in Eq. (18.5) using comoving-mesh adaptation
% aimed at equalizing the arclength across all collocation intervals in
% state space, as described in Sect. 18.3.1. Here, the distribution of mesh
% points is uniquely determined along the solution manifold by appending a
% family of zero functions to the continuation problem, corresponding to a
% mixed Euler method applied to a gradient vector field. In contrast to
% Examples 18.1 and 18.4, the comoving-mesh algorithm allows one to compute
% the solution family with the same desired tolerance as in Fig. 18.3, but
% with only 1/3 of the number of mesh points. Again, we need to restart
% twice with a finer mesh at terminal solution points. The solution
% obtained for epsilon=20 is shown in panel (a). Panel (b) shows the
% corresponding time profile and illustrates the adaptive time mesh. The
% estimated error plotted against variations in epsilon is shown in panel
% (c).

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run') ...
    && coco_exist('run3', 'run'))
  run demo_coll_v4
end

coco_use_recipes_toolbox po_v3 coll_v4

% Extract data
bd1  = coco_bd_read('run1');
bd2  = coco_bd_read('run2');
bd3  = coco_bd_read('run3');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-3.0833 3.0833 -1.0794*(1-1/6) 1.0794*(1-1/6)]*1.05)

lab  = coco_bd_labs(bd3, 'EP');   % Extract label
sol3 = po_read_solution('', 'run3', lab(end)); % Extract solution

idx = 1:5:50*5+1;
plot(sol3.x(idx,1), sol3.x(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(sol3.x(:,1), sol3.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel(b)
figure(2)
clf
hold on
grid on
box on

idx = 1:5:50*5+1;
plot(sol3.t(idx)/sol3.t(end), sol3.x(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(sol3.t/sol3.t(end), sol3.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol3.t(idx)/sol3.t(end), 0*sol3.x(idx,2)-0.95, 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

axis('tight')
lims = get(gca, 'YLim');
axis([-0.005 1.005 lims(1)-0.05*(lims(2)-lims(1)) ...
  lims(2)+0.05*(lims(2)-lims(1))])

hold off

% Plot data: panel (c)
figure(3)
clf
hold on
grid on
box on
axis([0 20 0 0.001])

p   = coco_bd_col(bd1, 'eps'); % Extract column data
err = coco_bd_col(bd1, 'po.seg.coll.err'); % Extract column data
plot(p, err, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

p   = coco_bd_col(bd2, 'eps'); % Extract column data
err = coco_bd_col(bd2, 'po.seg.coll.err'); % Extract column data
plot(p, err, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

p   = coco_bd_col(bd3, 'eps'); % Extract column data
err = coco_bd_col(bd3, 'po.seg.coll.err'); % Extract column data
plot(p, err, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

coco_use_recipes_toolbox

end
