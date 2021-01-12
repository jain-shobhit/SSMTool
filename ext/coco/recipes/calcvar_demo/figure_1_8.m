function figure_1_8
% Figure 1.8: Sample approximate extremal curves of the integral functional
% in Eq. (1.3) with action integrand given in Eq. (1.2) and satisfying the
% given boundary conditions at x=0 and x=1. The results shown here were
% obtained using a collocation method with N=10 and m=4 applied to a
% two-point boundary-value problem (a) as defined in Sect. 1.3.2, and to a
% quadrature approximation (b) as defined in Sect. 1.3.3. Each panel shows
% pairs of stationary curves for each given value of Y.

% Generate data
if ~(coco_exist('run_diff', 'run') && coco_exist('run_quad', 'run'))
  run demo_differential
  run demo_quadrature
end

coco_use_recipes_toolbox coll_v1 bvp_v1

% Extract data
bd   = coco_bd_read('run_diff');
labs = coco_bd_labs(bd, 'UZ');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-0.02 1.02 -0.02 1.02])

for lab=labs
  if lab==labs(2) || lab==labs(5)
    sol = bvp_read_solution('', 'run_diff', lab); % Extract solution
    plot(sol.t, sol.x(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
      'Color', [0.5 0.5 0.5], 'Marker', '.', 'MarkerSize', 15)
  else
    sol = bvp_read_solution('', 'run_diff', lab); % Extract solution
    plot(sol.t, sol.x(:,1), 'LineStyle', '-', 'LineWidth', 2, ...
      'Color', 'black', 'Marker', '.', 'MarkerSize', 15)
  end
end

hold off

% Extract data
bd   = coco_bd_read('run_quad');
labs = coco_bd_labs(bd, 'UZ');

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([-0.02 1.02 -0.08 4.08])

for lab=labs
  if lab==labs(2) || lab==labs(5)
    [calc chart] = coco_read_solution('calcvar', 'run_quad', lab); % Extract solution
    plot(calc.tbp, chart.x(calc.x_idx), 'LineStyle', '-', ...
      'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'Marker', '.', ...
      'MarkerSize', 15)
  else
    [calc chart] = coco_read_solution('calcvar', 'run_quad', lab); % Extract solution
    plot(calc.tbp, chart.x(calc.x_idx), 'LineStyle', '-', ...
      'LineWidth', 2, 'Color', 'black', 'Marker', '.', 'MarkerSize', 15)
  end
end

hold off

coco_use_recipes_toolbox

end
