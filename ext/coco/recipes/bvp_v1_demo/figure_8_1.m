function figure_8_1
% Figure 8.1: Solutions to the boundary-value problem given by the system
% of differential equations in Eqs. (8.6)-(8.7) and the given boundary
% conditions at x=0 and x=1. We start with the function y_1==1 shown in
% panel (a) as an initial guess. The initial correction converges to
% solution 1. A family of solutions obtained using numerical continuation
% is shown in panel (b). The labels correspond to the session output
% included in the text.

% Generate data
if coco_exist('catn', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_catenary

coco_use_recipes_toolbox coll_v1 bvp_v1

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0 1 0 3])

plot(t0, x0(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.3 0.3 0.3]) %#ok<NODEF>
sol = bvp_read_solution('', 'catn', 1); % Extract solution
plot(sol.t, sol.x(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 1 0 3])

bd   = coco_bd_read('catn'); % Extract bifurcation data
labs = coco_bd_labs(bd);     % Extract solution labels
for lab=labs
  sol = bvp_read_solution('', 'catn', lab); % Extract solution
  plot(sol.t, sol.x(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
end

hold off

coco_use_recipes_toolbox

end
