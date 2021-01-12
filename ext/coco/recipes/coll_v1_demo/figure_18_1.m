function figure_18_1
% Figure 18.1: Comparison of approximate periodic orbits of the dynamical
% system with vector field given by Eq. (18.5) at epsilon=20, obtained in
% Example 18.1 using continuation for different uniform meshes. For the
% default setting with N=10 and m=4, we obtain a poor approximation.
% Increasing the discretization parameters to N=150 and m=5 results in a
% solution curve that seems accetable, at least by visual inspection. A
% large subset of the additional mesh points contributes little to the
% improvement, however; only a small number is allocated along the
% problematic vertical part of the orbit.

% Generate data
if ~(coco_exist('pneta1', 'run') && coco_exist('pneta2', 'run'))
  run demo_pneta
end

coco_use_recipes_toolbox po_v1 coll_v1

% Extract data
bd  = coco_bd_read('pneta1');
lab = coco_bd_labs(bd, 'EP');
sol1 = po_read_solution('', 'pneta1', lab(end));
bd  = coco_bd_read('pneta2');
lab = coco_bd_labs(bd, 'EP');
sol2 = po_read_solution('', 'pneta2', lab(end));

% Plot data
figure(1)
clf
hold on
box on
grid on
axis([-3.0833 3.0833 -1.0794 1.0794]*1.05)

plot(sol2.x(:,1), sol2.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12)
plot(sol1.x(:,1), sol1.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

coco_use_recipes_toolbox

end
