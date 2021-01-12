function figure_19_2
% Figure 19.2: Continuation of periodic orbits of the dynamical system
% given by the vector field in Eq. (18.5) using the adaptive spectral
% toolbox, as implemented in Sects. 19.2 and 19.3. The rapid decay of
% Fourier coefficients illustrated in Fig. 19.1 enables adaptation by
% exchanging equations that either fix a avaiable at 0 or make it equal to
% a Fourier mode. In contrast to the results shown in Figs. 18.3 and 18.5,
% continuation here computes the complete family in a single run. Panel (a)
% shows sample orbits. The effects of adaptation are illustrated in panels
% (b) and (c). Panel (b) differs from the figure in the book. The latter
% appears to have been generated using an earlier version of the data.

% Note that machine-dependent round-off errors may result in different
% error estimates and discretization orders.

% Generate data
if ~(coco_exist('run1', 'run') && coco_exist('run2', 'run'))
  run demo_pneta
end

coco_use_recipes_toolbox dft_v1

% Extract data
[sol1 data] = dft_read_solution('', 'run1', 3);
sol2 = dft_read_solution('', 'run1', 6);
sol3 = dft_read_solution('', 'run1', 8);
dft = data.dft;
TOL1 = dft.TOLINC;
TOL2 = dft.TOLDEC;

bd1 = coco_bd_read('run1');
bd2 = coco_bd_read('run2');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
box on
axis([-3.0833 3.0833 -1.0794*(1-1/6) 1.0794*(1-1/6)]*1.05)

plot(sol3.x(:,1), sol3.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol2.x(:,1), sol2.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
plot(sol1.x(:,1), sol1.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 20 0 dft.TOL])

plot([0 20], [TOL1 TOL1], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', [0.3 0.3 0.3])
plot([0 20], [TOL2 TOL2], 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', [0.3 0.3 0.3])

p   = coco_bd_col(bd1, 'eps');
err = coco_bd_col(bd1, 'dft.err');
plot(p, err, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)
p   = coco_bd_col(bd2, 'eps');
err = coco_bd_col(bd2, 'dft.err');
plot(p, err, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panel (c)
figure(3)
clf
hold on
grid on
box on
axis([0 20 0 300])

nm = coco_bd_col(bd1, 'dft.NMOD');
p  = coco_bd_col(bd1, 'eps');
plot(p, nm, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)

nm = coco_bd_col(bd2, 'dft.NMOD');
p  = coco_bd_col(bd2, 'eps');
plot(p, nm, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)

hold off

coco_use_recipes_toolbox

end
