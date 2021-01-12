function figure_7_1
% Figure 7.1: A technique for constructin initial solutions to nonlinear
% boundary-value problems is continuation in the interval length T. Here,
% one starts with a short trajectory segment for small T<<1 and continues
% in T until reaching the desired value. This process is called growing an
% initial orbit and is illustrated in Sect. 7.3.1 in the context of the
% problem from calculus of variations, introduced in Chap. 1. Panel (a)
% shows an initial guess obtained with an Euler step at t=0 with step size
% h=0.04 and compares this approximation with the exact solution. The grown
% solution is shown in panel (b). From this solution one can, finally,
% start a continuation in Y. The labeled solution curves in (c) correspond
% to the labels in the session output included in the text.

% Generate data
if coco_exist('coll1', 'run') && coco_exist('coll2', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_catenary

coco_use_recipes_toolbox coll_v1

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([0 1 0 2])

t = linspace(0,1,1000);
x = cosh(t);

plot(t, x, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
plot(t0, x0(:,1), 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12) %#ok<NODEF>

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 1 0 2])

sol = coll_read_solution('', 'coll1', 5); % Extract solution
plot(sol.t, sol.x(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (c)
figure(3)
clf
hold on
grid on
box on
axis([0 1 0 3])

bd = coco_bd_read('coll2'); % Extract bifurcation data
labs = coco_bd_labs(bd);    % Extract solution labels
for lab=labs
  sol = coll_read_solution('', 'coll2', lab); %Extract solution
  plot(sol.t, sol.x(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

end

hold off

coco_use_recipes_toolbox

end
