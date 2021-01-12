function figure_8_6
% Figure 8.6: Continuation of the periodic orbit with high period,
% illustrated in Fig. 8.5, while keeping the period constant, resulting in
% an approximation to a homoclinic bifurcation curve (a). Each point on
% this curve corresponds to a terminal point along a family of periodic
% orbits emenating from a Hopf bifurcation under variations in p_2. Panel
% (b) shows selected members of the family of high-period orbits. The
% labels correspond to the session output included in the text.

% Generate data
if coco_exist('po1', 'run') && coco_exist('po2', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_marsden

coco_use_recipes_toolbox coll_v1 po_v1

% Extract data
bd   = coco_bd_read('po2');
p1 = coco_bd_col(bd, 'p1');
p2 = coco_bd_col(bd, 'p2');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-0.08 0.0 1 21])

plot(p1, p2, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([-0.7 0.4 -0.3 1.3])

labs = [1 3 6 9 12];
for lab=labs;
  sol = po_read_solution('', 'po2', lab); % Extract solution
  plot(sol.x(:,2), sol.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
end

hold off

coco_use_recipes_toolbox

end
