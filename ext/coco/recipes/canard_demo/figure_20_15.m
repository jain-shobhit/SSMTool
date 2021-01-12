function figure_20_15
% Figure 20.15: Selected approximate periodic orbits obtained during
% continuation along the canard family of the dynamical system given by the
% vector field in Eq. (20.45). Both the adaptive spectral method (a) from
% Sect. 19.3 with NMAX equal to 100 and the nonadaptive collocation method
% (b) from Sect. 18.1.2 with NTST and NCOL equal to 100 and 4,
% respectively, fail early on the canard family, given the desired error
% tolerance of 10^(-2).

% Generate data
if ~(coco_exist('1', 'run') && coco_exist('2', 'run'))
  run demo_canard
end

coco_use_recipes_toolbox dft_v1

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-1 0.8 -3 3])

bd   = coco_bd_read('1');
labs = coco_bd_labs(bd, 'UZ');
lab  = coco_bd_labs(bd, 'MXCL');
sol  = dft_read_solution('', '1', lab); % Extract solution

plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.4 0.4 0.4], 'Marker', '.', 'MarkerSize', 12)
for lab=labs
  sol = dft_read_solution('', '1', lab); % Extract solution
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
end

hold off

coco_use_recipes_toolbox

coco_use_recipes_toolbox po_v1 coll_v3

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([-1 0.8 -3 3])

bd   = coco_bd_read('2');
labs = coco_bd_labs(bd, 'UZ');
lab = coco_bd_labs(bd, 'MXCL');
sol = po_read_solution('', '2', lab); % Extract solution

plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', [0.4 0.4 0.4], 'Marker', '.', 'MarkerSize', 12)
for lab=labs
  sol = po_read_solution('', '2', lab); % Extract solution
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
end

hold off

coco_use_recipes_toolbox

end
