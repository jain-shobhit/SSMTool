function figure_8_4
% Figure 8.4: The last periodic orbit 6 obtained in the continuation run
% shown in Fig 8.3. The left panel shows the phase plot and the right panel
% the time profile. The time profile shows a phase of slow dynamics,
% indicating existence of a nearby equilibrium point, followed by a fast
% excursion. To confirm our hypothesis, we use this orbit and extent the
% part that seems to be close to an equilibrium; see Fig. 8.5.

% Generate data
if coco_exist('po2', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_marsden

coco_use_recipes_toolbox coll_v1 po_v1

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-0.5 0.1 -0.1 0.8])

sol = po_read_solution('', 'po1', 6); % Extract solution
plot(sol.x(:,2), sol.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([-0.5 20 -0.1 0.8])

plot(sol.t, sol.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

coco_use_recipes_toolbox

end
