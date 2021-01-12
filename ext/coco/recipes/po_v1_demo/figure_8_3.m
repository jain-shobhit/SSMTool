function figure_8_3
% Figure 8.3: We compute a family of periodic orbits emanating from a Hopf
% bifurcation point of the dynamical system given by the vector field in
% Eq. (8.24). We obtain an initial solution guess from normal form
% analysis, shown in gray in panel (a). The initial correction step
% coverges to orbit 1. Panel (b) shows the family of periodic orbits of
% increasing amplitudes that seem to approach a homoclinic orbit, indicated
% by the corner that develops in the top left part of the plot and
% allocates many mesh points due to slow dynamics. The labels correspond to
% the session output included in the text.

% Generate data
if coco_exist('po1', 'run')
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
axis(0.01*[-1.25 1.25 -1.25 1.25])

t0 = (0:2*pi/100:2*pi)';
x0 = 0.01*(cos(t0)*[1 0 -1]-sin(t0)*[0 1 0]);

plot(x0(:,2), x0(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.5 0.5 0.5])
sol = po_read_solution('', 'po1', 1); % Extract solution
plot(sol.x(:,2), sol.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)

figure(2)
clf
hold on
grid on
box on
axis([-0.5 0.1 -0.1 0.8])

labs = 1:6;
for lab=labs
  sol = po_read_solution('', 'po1', lab); % Extract solution
  plot(sol.x(:,2), sol.x(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
end

hold off

coco_use_recipes_toolbox

end
