function figure_8_7
% Figure 8.7: Continuation of periodic orbits of the Lienard system in
% Example 8.5. Both continuation runs start at the same initial solution,
% marked with label 1. The initial point y(0) on each orbit is emphasized.
% With a fixed phase condition, all initial points must lie on a straight
% line (a), and we are unable to compute the full family due to a tangency.
% Updating the phase condition in each continuation step results in the
% full family; the inintial points now lie on a curve passing through the
% Hopf bifurcation point at the origin within numerical accuracy (b).

% Generate data
if ~coco_exist('run_fixed', 'run') || ~coco_exist('run_moving', 'run')
  run demo_lienard
end

coco_use_recipes_toolbox coll_v1 bvp_v2

% Plot data: panel (a)
bd = coco_bd_read('run_fixed'); % Extract bifurcation data
labs = coco_bd_labs(bd);        % Extract labels
x0 = [];

figure(1)
clf
hold on
grid on
box on
axis([-1.75 1.75 -1.4 1.23])

for lab=labs
  sol = bvp_read_solution('', 'run_fixed', lab);
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.7 0.7 0.7], 'Marker', '.', 'MarkerSize', 12)
  x0 = [x0 ; sol.x(1,:)]; %#ok<AGROW>
end
plot(x0(:,1), x0(:,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15)

hold off

% Plot data: panel (b)
bd = coco_bd_read('run_moving'); % Extract bifurcation data
labs = coco_bd_labs(bd);         % Extract labels
x0 = [];

figure(2)
clf
hold on
grid on
box on
axis([-1.75 1.75 -1.4 1.23])

for lab=labs
  sol = bvp_read_solution('', 'run_moving', lab);
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.7 0.7 0.7], 'Marker', '.', 'MarkerSize', 12)
  x0 = [x0 ; sol.x(1,:)]; %#ok<AGROW>
end
plot(x0(:,1), x0(:,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15)

hold off

coco_use_recipes_toolbox

end
