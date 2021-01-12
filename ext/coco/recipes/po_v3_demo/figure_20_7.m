function figure_20_7
% Figure 20.7: State-space representations of near-homoclinic approximate
% periodic orbits of the dynamical system given by the vector field in Eq.
% (20.44), obtained using continuation with the nonadaptive 'coll' toolbox
% from Sect. 18.1.2 with NTST equal to 50 and TOL equal to 10^(-4). Panels
% (a) and (b) show the terminal low-period solution (gray circles) with 50
% mesh intervals found using continuation with a 1-dimensional atlas
% algorithm starting at a Hopf bifurcation point, as well as the corrected
% high-period solution (black curve and end points of mesh intervals
% denoted by x's) found using continuation with a 0-dimensional atlas
% algorithm applied to a reconstructed periodic orbit with 1,250 mesh
% intervals. The estimated error is 4.2532*10^(-5) for the low-period orbit
% and 4.2531*10^(-5) for the high-period orbit. These results are compared
% in Figs. 20.9 and 20.11 with computations using moving-mesh adaptation.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run'))
  run demo_po_v3
end

coco_use_recipes_toolbox po_v1 coll_v3

% Extract data
sol0 = po_read_solution('', 'run1', 6);
sol1 = po_read_solution('', 'run2', 1);

idx1 = 1:4:50*4+1;
idx2 = 1:4:1250*4+1;

% Plot data: panel (a)
figure(1)
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

% Plot data: panel (b)
figure(2)
clf
hold on
box on
grid on
axis([-0.5 0.1 -0.1 0.8])

plot(sol0.x(idx1,2), sol0.x(idx1,3), 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', [0.4 0.4 0.4], 'Marker', 'o', ...
  'MarkerSize', 8, 'MarkerFaceColor', 'white')
plot(sol1.x(:,2), sol1.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol1.x(idx2,2), sol1.x(idx2,3), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'x', 'MarkerSize', 9.5)

hold off

coco_use_recipes_toolbox

end
