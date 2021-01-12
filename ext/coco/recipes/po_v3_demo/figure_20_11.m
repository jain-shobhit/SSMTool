function figure_20_11
% Figure 20.11: Mesh adaptation for an approximate near-homoclinic periodic
% orbit of the dynamical system given by the vector field in Eq. (20.44),
% obtained using the collocation method with a moving mesh of variable
% discretization order from Sect. 20.2.2 with TOL equal to the default
% tolerance of 10^(-4). Panels (a) and (b) show variations in the
% discretization error and discretization order during 100 iterations of a
% 0-dimensional remesh-correct cycle applied to the reconstructed initial
% solution guess obtained by a 5,000-fold increase in the period T and a
% 6-fold increase in discretization order N. Again, we observe that the
% iteration settles onto an acceptable mesh after a somewhat longer
% transient phase. Although the approximation error drops much faster than
% in Fig. 20.9(a), the dynamics is dominated by the slow reduction of the
% number of mesh intervals. After 40 steps, the solution settles onto a
% mesh of size N=79. A comparison of the resulting solutions is shown in
% Fig. 20.12. Fig 20.13 shows the adapted mesh and time profile of the
% high-period orbit.

% Note that machine-dependent round-off errors may result in differences in
% the discretization errors and orders.

% Generate data
if ~(coco_exist('run1b', 'run') && coco_exist('run2b', 'run'))
  run demo_po_v3
end

coco_use_recipes_toolbox po_v3 coll_v6

% Extract data
bd  = coco_bd_read('run2b');
pt  = coco_bd_col(bd,'PT');
err = log10(coco_bd_col(bd,'po.seg.coll.err'));
N   = coco_bd_col(bd,'po.seg.coll.NTST');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0 100 -7 -2])

plot(pt, err, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)
plot(pt, 0*err+log10(5.0e-5), 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', [0.3 0.3 0.3])
plot(pt, 0*err+log10(1.0e-5), 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', [0.3 0.3 0.3])

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 100 60 180])

plot(pt, N, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)

hold off

coco_use_recipes_toolbox

end
