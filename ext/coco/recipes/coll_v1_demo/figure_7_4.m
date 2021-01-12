function figure_7_4
% Figure 7.4: Computation of an initial approximation of a heteroclinic
% connection in the dynamical system with vector field given in Eq. (7.39),
% following the methodology described in Sect. 7.3.3; cf. Fig. 7.3. Here,
% we combine two instances of the 'coll' toolbox with two instances of the
% 'alg' toolbox, which are used to solve for the equilibria and an
% invariant eigenspace. A family of connecting orbits is shown in panel
% (f). The labels correspond to the session output of the last run included
% in the text.

% Generate data
if coco_exist('doedel1', 'run') && coco_exist('doedel2', 'run') ...
    && coco_exist('doedel3', 'run') && coco_exist('doedel4', 'run') ...
    && coco_exist('doedel5', 'run') && coco_exist('doedel6', 'run')
  coco = @(prob, run, varargin) coco_bd_read(run); %#ok<NASGU>
end
run demo_doedel

coco_use_recipes_toolbox coll_v1 alg_v6

compute_extra_data(segs, algs, eps0, th0, vec0, lam0);

% Plot data: panel (a)
N     = 30;
alims = [-1.2 1.5 -2 2.2];
x     = linspace(alims(1),alims(2),N);
y     = linspace(alims(3),alims(4),N);
[X Y] = meshgrid(x,y);
fxy   = doedel([X(:) Y(:)]', repmat([1;1], [1 numel(X)]));
fx    = reshape(fxy(1,:), size(X));
fy    = reshape(fxy(2,:), size(X));

figure(1)
clf
hold on
grid on
box on
axis([-1.2 1.5 -2 2.2])

quiver(X, Y, fx, fy, 2, 'LineStyle', '-', 'LineWidth',   1, ...
  'Color', [0.3 0.3 0.3], 'ShowArrowHead', 'on',  'MaxHeadSize', 0.2)
plot_skeleton('LineStyle', '-', 'LineWidth', 2, 'Color', 'black')
plot(-1, 1, 'LineStyle', '-', 'LineWidth', sqrt(2), 'Color', 'black', ...
  'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'white')

hold off

% Plot data: panels (b)-(e)
labs = [1, 3, 5, 2];
runs = {'doedel1', 'doedel1', 'doedel2', 'doedel5'};

for i=1:4
  figure(i+1)
  clf
  hold on
  grid on
  box on
  axis([-1.2 1.5 -2 2.2])
  
  plot_skeleton('LineStyle', '-', 'LineWidth', 2, 'Color', [0.7 0.7 0.7])
  
  sol = coll_read_solution('doedel1', runs(i), labs(i)); % Extract solution
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
  sol = coll_read_solution('doedel2', runs(i), labs(i)); % Extract solution
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12)
  
  hold off
end

% Plot data: panel (f)
figure(6)
clf
hold on
grid on
box on
axis([-1.2 1.5 -2 2.2])

bd = coco_bd_read('doedel6'); % Extract bifurcation data
labs = coco_bd_labs(bd);      % Extract solution labels
labs = labs(labs<=4 | labs==6 | labs==7 | labs==10);
for lab=labs
  sol = coll_read_solution('doedel1', 'doedel6', lab); % Extract solution
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
  sol = coll_read_solution('doedel2', 'doedel6', lab); % Extract solution
  plot(sol.x(:,1), sol.x(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
end

hold off

coco_use_recipes_toolbox

end

function compute_extra_data(segs, algs, eps0, th0, vec0, lam0)

% high precision runs for skeleton
if ~coco_exist('doedel1a', 'run')
  prob = coco_prob();
  prob = coco_set(prob, 'coll', 'NTST', 100);
  prob = doedel_isol2het(prob, segs, algs, eps0, th0, vec0, lam0);
  coco(prob, 'doedel1a', [], 1, 'y12e', [0 0.99]);
end
if ~coco_exist('doedel2a', 'run')
  prob = coco_prob();
  prob = coco_set(prob, 'cont', 'ItMX', 100);
  prob = doedel_sol2het(prob, 'doedel1a', 3);
  coco(prob, 'doedel2a', [], 1, 'y22e', [-0.995 0]);
end
if ~coco_exist('doedel3a', 'run')
  prob = doedel_sol2het(coco_prob(), 'doedel2a', 8);
  coco(prob, 'doedel3a', [], 1, 'gap', [-2 0]);
end
if ~coco_exist('doedel4a', 'run')
  prob = doedel_sol2het(coco_prob(), 'doedel3a', 6);
  coco(prob, 'doedel4a', [], 1, 'eps1', [1e-3 eps0(1)]);
end
if ~coco_exist('doedel5a', 'run')
  prob = doedel_sol2het(coco_prob(), 'doedel4a', 3);
  coco(prob, 'doedel5a', [], 1, 'eps2', [1e-3 eps0(2)]);
end

end

function plot_skeleton(varargin)

v1 = [-3/sqrt(10); 1/sqrt(10)];
[t x] = ode45(@(t,x) -doedel(x,[1;1]), [0 4.2], [1;-1]-1.0e-4*v1); %#ok<ASGLU>
plot(x(:,1), x(:,2), varargin{:})
sol = coll_read_solution('doedel1', 'doedel5a', 3);
plot([-1 ; sol.x(:,1)], [1 ; sol.x(:,2)], varargin{:});
sol = coll_read_solution('doedel2', 'doedel5a', 3);
plot([sol.x(:,1);x(1,1)], [sol.x(:,2);x(1,2)], varargin{:});
plot([1 1], [-3 3], varargin{:});

end
