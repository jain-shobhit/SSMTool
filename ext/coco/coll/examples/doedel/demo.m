%% Continuation of heteroclinic connections from Section 7.3.3 of Recipes for Continuation
%
% Computation of an initial approximation of a heteroclinic connection in
% the dynamical system with vector field given in Eq. (7.39) in Recipes for
% Continuation, following the methodology described in Sect. 7.3.3; cf.
% Fig. 7.3. Here, we combine two instances of the 'coll' toolbox with two
% instances of the 'alg' toolbox, which are used to solve for the
% equilibria and an invariant eigenspace. 

% A family of connecting orbits is shown in panel (f). The labels
% correspond to the session output of the last run included in the text.

%% Encoding

% The continuation problem encoded below includes seven monitor functions,
% and the corresponding inactive continuation parameters 'gap', 'p1', 'p2',
% 'epts1', 'epts2', 'y12e', and 'y22e'. Its dimensional deficit equals 0.
% The calls to the coco entry-point function indicate a desired manifold
% dimension of 1. To this end, one continuation parameters is released and
% allowed to vary during each continuation run.

p0 = [1; 1];
dev0 = [0.05; 0.05];
th0 = -pi/2;
eqs10 = [-1; 1];
eqs20 = [1; -1];
vec0 = [-3/sqrt(10); 1/sqrt(10)];
lam0 = -2;
segs(1).t0 = [0; 1];
x0         = eqs10+dev0(1)*[cos(th0); sin(th0)];
segs(1).x0 = [x0  x0+doedel(x0, p0)]';
segs(1).p0 = p0;
segs(2).t0 = [0; 1];
x0         = eqs20+dev0(2)*vec0;
segs(2).x0 = [x0-doedel(x0, p0) x0]';
segs(2).p0 = p0;
epts(1).x0 = eqs10;
epts(1).p0 = p0;
epts(2).x0 = eqs20;
epts(2).p0 = p0;

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 3);
prob = coco_set(prob, 'corr', 'TOL', 1e-7);

fprintf('\n Run=''%s'': Continue trajectory segments until y12e=0\n', ...
  'doedel1');

pprob = doedel_isol2het(prob, segs, epts, dev0, th0, vec0, lam0);
coco(pprob, 'doedel1', [], 1, {'y12e'  'doedel1.coll.err_TF'}, [0 0.99]);

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s'' until y22e=0\n', ...
  'doedel2', 3, 'doedel1');

pprob = doedel_sol2het(prob, 'doedel1', 3);
coco(pprob, 'doedel2', [], 1, {'y22e' 'doedel1.coll.err_TF'}, [-0.995 0]);

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s'' until gap=0\n', ...
  'doedel3', 6, 'doedel2');

pprob = doedel_sol2het(prob, 'doedel2', 6);
coco(pprob, 'doedel3', [], 1, {'gap'  'doedel1.coll.err_TF'}, [-2 0]);

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s'' until dev1=1e-3\n', ...
  'doedel4', 5, 'doedel3');

pprob = doedel_sol2het(prob, 'doedel3', 5);
coco(pprob, 'doedel4', [], 1, {'dev1'  'doedel1.coll.err_TF'}, [1e-3 dev0(1)]);

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s'' until dev2=1e-3\n', ...
  'doedel5', 3, 'doedel4');

pprob = doedel_sol2het(prob, 'doedel4', 3);
coco(pprob, 'doedel5', [], 1, {'dev2'  'doedel1.coll.err_TF'}, [1e-3 dev0(2)]);

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s''\n', ...
  'doedel6', 2, 'doedel5');

pprob = doedel_sol2het(prob, 'doedel5', 2);
coco(pprob, 'doedel6', [], 1, {'p2'  'doedel1.coll.err_TF'}, [0.5 8]);

%% Graphical representation of stored solutions


% Plot data: panel (a)
N     = 30;
alims = [-1.2 1.5 -2 2.2];
x     = linspace(alims(1),alims(2),N);
y     = linspace(alims(3),alims(4),N);
[X, Y] = meshgrid(x,y);
fxy   = doedel([X(:) Y(:)]', repmat([1;1], [1 numel(X)]));
fx    = reshape(fxy(1,:), size(X));
fy    = reshape(fxy(2,:), size(X));

figure(1); clf; hold on; grid on; box on; axis([-1.2 1.5 -2 2.2])

quiver(X, Y, fx, fy, 2, 'LineStyle', '-', 'LineWidth',   1, ...
  'Color', [0.3 0.3 0.3], 'ShowArrowHead', 'on',  'MaxHeadSize', 0.2)
plot_skeleton('LineStyle', '-', 'LineWidth', 2, 'Color', 'black')
plot(-1, 1, 'LineStyle', '-', 'LineWidth', sqrt(2), 'Color', 'black', ...
  'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'white')

hold off

% Plot data: panels (b)-(e)
labs = [1, 3, 6, 2];
runs = {'doedel1', 'doedel1', 'doedel2', 'doedel5'};

for i=1:4
  figure(i+1); clf; hold on; grid on; box on  
  plot_skeleton('LineStyle', '-', 'LineWidth', 2, 'Color', [0.7 0.7 0.7])
  coco_plot_sol(runs{i}, labs(i), 'doedel', 1:2, 'x', 'x')
  hold off; axis([-1.2 1.5 -2 2.2])
end

% Plot data: panel (f)
figure(6); clf; hold on; grid on; box on

bd = coco_bd_read('doedel6'); % Extract bifurcation data
labs = coco_bd_labs(bd);      % Extract solution labels
labs = labs(labs<=4 | labs==6 | labs==7 | labs==10);
for lab=labs
  coco_plot_sol('doedel6', lab, 'doedel', 1:2, 'x', 'x')
  sol = coll_read_solution('doedel1', 'doedel6', lab); % Extract solution
  v1 = [-(-2 + sol.p(2))/sol.p(1), 1];
  v1 = v1/norm(v1);
  v2 = [0, 1];
  plot([-1,-1+0.01*v1(1)], [sol.p(1)/sol.p(2), sol.p(1)/sol.p(2)+0.01*v1(2)])
  plot([-1,-1+0.01*v2(1)],[sol.p(1)/sol.p(2), sol.p(1)/sol.p(2)+0.01*v2(2)])
  axis([-1 1 -2 2.2])
  pause(0.1)
  
end
hold off
