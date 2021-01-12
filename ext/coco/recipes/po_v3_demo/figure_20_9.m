function figure_20_9
% Figure 20.9: Mesh adaptation and an approximate near-homoclinic periodic
% orbit of the dynamical system given by the vector field in Eq. (20.44),
% obtained using the collocation method with a moving mesh of fixed order
% from Sect. 20.2.1 with TOL equal to 10^(-4). Panel (a) shows variations
% in the discretization error during 100 iterations of a 0-dimensional
% remesh-correct cycle applied to the reconstructed initial solution guess
% obtained by a 5,000-fold increase in the period T and a 6-fold increase
% in discretization order N. In panels (b) and (c), we compare the
% corrected high-period orbit (black) with the low-period orbit (gray
% circles) obtained using continuation from the Hopf bifurcation point. The
% number of mesh points along the excursion from the equilibrium is nearly
% identical for both orbits, although the mesh points move somewhat. The
% estimated error is 2.2619*10^(-5) for the low-period orbit and
% 4.2620*10^(-5) for the high-period orbit. The adapted mesh and time
% profile of the high-period orbit are shown in Fig. 20.10. Figs.
% 20.11-20.13 repeat the analysis with a moving mesh with variable
% discretization order.

% Note that TOL is set to 10^(-1) during continuation of the high-period
% orbits. Machine-dependent round-off errors may affect the discretization
% errors.

% Generate data
if ~(coco_exist('run1a', 'run') && coco_exist('run2a', 'run'))
  run demo_po_v3
end

coco_use_recipes_toolbox po_v3 coll_v5

% Extract data
bd = coco_bd_read('run2a');

sol0 = po_read_solution('', 'run1a', 6);
sol1 = po_read_solution('', 'run2a', 2);
                                         
idx1 = 1:4:20*4+1;
idx2 = 1:4:120*4+1;

% Plot data: panel (a)
figure(1)
clf
hold on
box on
grid on

pt  = coco_bd_col(bd,'PT');
err = log10(coco_bd_col(bd,'po.seg.coll.err'));
plot(pt, err, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)
plot(pt, 0*err-4, 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', [0.3 0.3 0.3])

axis('tight')  
ylims = get(gca, 'YLim');
axis([0 100 ylims(1)-0.05*(ylims(2)-ylims(1)) ...
  ylims(2)+0.05*(ylims(2)-ylims(1))])

hold off

% Plot data: panel (b)
figure(2)
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
plot(sol1.x(idx2,1), sol1.x(idx2,2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'x', 'MarkerSize', 9.5)

hold off

% Plot data: panel (c)
figure(3)
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
