function figure_20_14
% Figure 20.14: A canard explosion in the Van der Pol dynamical system
% given by the vector field in Eq. (20.45). Here, a dramatic increase of
% amplitude (represented by the norm ||U||=||(u,mu_J)||) along a family of
% periodic orbits over an exceedingly small variation in a system parameter
% shows up as a virtually vertical branch in a bifurcation diagram (a).
% Panels (b) to (i) show members of the canard family. The orbits along the
% family were selected according to their period, indicated in the caption.
% We observe evidence for a fold point with respect to period between
% orbits 5 and 6.

% Generate data
if ~(coco_exist('4', 'run'))
  run demo_canard
end

coco_use_recipes_toolbox po_v3 coll_v6

% Extract data
bd   = coco_bd_read('4');
labs = coco_bd_labs(bd, 'UZ');
a    = coco_bd_col(bd, 'a');
y    = coco_bd_col(bd, '||U||');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-1.01 -0.9 85 700])

plot(a, y, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black')

hold off

% Plot data: panels (b)-(i)
for i=1:8
  figure(i+1)
  clf
  hold on
  grid on
  box on
  axis([-1 0.8 -3 3])

  sol = po_read_solution('', '4', labs(i)); % Extract solution
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black')
  
  hold off
end

coco_use_recipes_toolbox

end
