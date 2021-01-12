function figure_8_5
% Figure 8.5: Starting with orbit 6 from Fig. 8.4, we insert a long segment
% of constant dynamics and rescale the period such that the shape of the
% orbit in phase space should be unchanged if there exists a nearby
% homoclinic orbit. The extended time profile after the initial correction
% step is shown in panel (a). We clearly observe an elongated phase of
% near-constant dynamics. We overlay this new solution 1 (black dot) on top
% of the previous orbit 6 (gray circle) in panel (b). The phase plots,
% including the distribution of mesh points, are virtually identical, which
% supports the assumption that a nearby homoclinic orbit exists. We
% continue a family of high-period orbit in Fig. 8.6.

% Generate data
if coco_exist('po1', 'run') && coco_exist('po2', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_marsden

coco_use_recipes_toolbox coll_v1 po_v1

% Extract data
sol0 = po_read_solution('', 'po1', 6);
sol1 = po_read_solution('', 'po2', 1);

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-10 500 -0.1 0.8])

plot(sol1.t, sol1.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([-0.5 0.1 -0.1 0.8])

plot(sol0.x(:,2), sol0.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.4 0.4 0.4], 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerFaceColor', 'white')
plot(sol1.x(:,2), sol1.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15)

hold off

coco_use_recipes_toolbox

end
