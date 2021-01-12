function figure_20_17
% Figure 20.17: Discretization of the canard orbit in Fig. 20.14(g)
% obtained using continuation with the moving-mesh method from Sect. 20.2.1
% on a mesh with NTST equal to 70, NCOL equal to 4, and TOL equal to
% 10^(-4). The fat dots in (a) and the open circles in (b) mark the end
% points of mesh intervals. The estimated error of this solution is
% 5.1787*10^(-6).

% Note that machine-dependent round-off errors may affect the
% discretization errors.

% Generate data
if ~(coco_exist('3', 'run') && coco_exist('4', 'run') ...
     && coco_exist('5', 'run'))
   run demo_canard
end

coco_use_recipes_toolbox po_v3 coll_v5

% Extract data
bd   = coco_bd_read('3');
labs = coco_bd_labs(bd, 'UZ');
lab  = labs(6);
sol  = po_read_solution('', '3', lab);

N   = (numel(sol.t)-1)/4;
idx = 1:4:N*4+1;

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-1 0.8 -3 3])

plot(sol.x(idx,1), sol.x(idx,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')
plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
box on
grid on
axis([0 1 -2.5 2])

tt  = sol.t/sol.t(end);
plot(tt, sol.x, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(tt(idx), 0*sol.x(idx,:)-2.3, 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white')

hold off

coco_use_recipes_toolbox

end
