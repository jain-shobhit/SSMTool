function figure_20_5
% Figure 20.5: Continuation of periodic orbits of the dynamical system
% given by the vector field in Eq. (18.5) using an equidistributed moving
% mesh without adaptive changes to the discretization order, as described
% in Sect. 20.2.1. The algorithm allows one to compute the solution family
% with the same desired tolerance as in Fig. 18.3 with 1/9 of the number of
% mesh points. Panel (a) shows sample orbits obtained during continuation,
% and panels (b) and (c) illustrate the effect of adaptation. A refined
% strategy with varying discretization order is illustrated in Fig. 20.6.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run'))
  run demo_coll_v5
end

coco_use_recipes_toolbox po_v3 coll_v5

% Extract data
bd1  = coco_bd_read('run1');
bd2  = coco_bd_read('run2');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-3.0833 3.0833 -1.0794*(1-1/6) 1.0794*(1-1/6)]*1.05)

lab = coco_bd_labs(bd1, 'MXCL'); % Extract label
sol1 = po_read_solution('', 'run1', lab); % Extract solution
lab = coco_bd_labs(bd2, 'EP');   % Extract label
sol2 = po_read_solution('', 'run2', lab(end)); % Extract solution

idx = 1:4:10*4+1;
plot(sol1.x(idx,1), sol1.x(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(sol1.x(:,1), sol1.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

idx = 1:4:20*4+1;
plot(sol2.x(idx,1), sol2.x(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(sol2.x(:,1), sol2.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
box on
grid on

idx = 1:4:20*4+1;
plot(sol2.t(idx)/sol2.t(end), sol2.x(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(sol2.t/sol2.t(end), sol2.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol2.t(idx)/sol2.t(end), 0*sol2.x(idx,2)-0.95, 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

axis('tight')
xlims = get(gca, 'XLim');
ylims = get(gca, 'YLim');
axis([xlims(1)-0.005*(xlims(2)-xlims(1)) ...
  xlims(2)+0.005*(xlims(2)-xlims(1)) ...
  ylims(1)-0.05*(ylims(2)-ylims(1)) ylims(2)+0.05*(ylims(2)-ylims(1))])

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

p   = coco_bd_col(bd2, 'eps');
err = coco_bd_col(bd2, 'po.seg.coll.err');
plot(p, err, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

coco_use_recipes_toolbox

end
