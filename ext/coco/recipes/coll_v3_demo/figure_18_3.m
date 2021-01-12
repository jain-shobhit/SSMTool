function figure_18_3
% Figure 18.3: Continuation of periodic orbits of the dynamical system
% given by the vector field in Eq. (18.5) using pointwise adaptive changes
% to the discretization order, as described in Sect. 18.1.2. Here, a
% terminal event is associated with a monitor function that computes an
% estimate of the discretization error. Continuation stops automatically
% whenever the estimated error exceeds a predefined limit. As shown in
% Example 18.4, this approach allows one to compute an approximate solution
% family within a desired tolerance. Here, we need to restart twice with a
% finer mesh at terminal solution points, represented by the corresponding
% approximate periodic orbits in panel (a) together with the solution for
% epsilon=20. The reason for the failure to remain within the desired
% tolerance with low discretization order is explained by the shart fronts
% evident in the corresponding time profiles in panel (b). Panel (c) shows
% the estimated error plotted against variations in epsilon.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run') ...
    && coco_exist('run3', 'run'))
  run demo_coll_v3
end

coco_use_recipes_toolbox po_v1 coll_v3

% Extract data
bd1  = coco_bd_read('run1');
bd2  = coco_bd_read('run2');
bd3  = coco_bd_read('run3');

% Plot data: panel (a)
figure(1)
clf
hold on
box on
grid on
axis([-3.0833 3.0833 -1.0794*(1-1/6) 1.0794*(1-1/6)]*1.05)

lab  = coco_bd_labs(bd1, 'MXCL'); % Extract label
sol1 = po_read_solution('', 'run1', lab); % Extract solution
lab  = coco_bd_labs(bd2, 'MXCL'); % Extract label
sol2 = po_read_solution('', 'run2', lab); % Extract solution
lab  = coco_bd_labs(bd3, 'EP');   % Extract label
sol3 = po_read_solution('', 'run3', lab(end)); % Extract solution

plot(sol3.x(:,1), sol3.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol2.x(:,1), sol2.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol1.x(:,1), sol1.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on

plot(sol3.t/sol3.t(end), sol3.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol2.t/sol2.t(end), sol2.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol1.t/sol1.t(end), sol1.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

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

p   = coco_bd_col(bd1, 'eps');
err = coco_bd_col(bd1, 'po.seg.coll.err');
plot(p, err, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

p   = coco_bd_col(bd2, 'eps');
err = coco_bd_col(bd2, 'po.seg.coll.err');
plot(p, err, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

p   = coco_bd_col(bd3, 'eps');
err = coco_bd_col(bd3, 'po.seg.coll.err');
plot(p, err, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

coco_use_recipes_toolbox

end
