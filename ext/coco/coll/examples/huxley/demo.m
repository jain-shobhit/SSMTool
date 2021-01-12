%% Continuation of heteroclinic connections from Section 7.3.2 of Recipes for Continuation
%
% In contrast to the example in Section 7.3.1 of Recipes for Continuation,
% the construction of an initial solution for the continuation of a
% heteroclinic orbit in the Huxley dynamical system given by the vector
% field in Eq. (7.26) consists of a sequence of steps. We start with two
% short orbit segments close to each saddle, initial guesses to which can
% be obtained by Euler steps in the respective invariant eigenspace. The
% two segments after initial correction are shown in Figure(2). In the
% subsequent step we grow both segments until they terminate on the
% hyperplane y_1=0.5, as shown in Figures 3 and 4. At this point, the two
% end points are separated by a small gap, the so-called Lin gap, which we
% close in the next step, as shown in Figure 5. In the last step, depicted
% in Figures 6 and 7, we grow the two segments toward the saddle equilibria
% and obtain an approximation to a connecting orbit, which can be used as a
% starting point for subsequent continuation runs.

% A skeleton of the dynamics of the Huxley dynamical system given by the
% vector field in Eq. (7.26) and its direction field for p_1=1/2 and p_2=0
% is shown in Figure 1. The system has three equilibria: a center at
% (1/2,0) and two saddles that are connected in a heteroclinic cycle.

%% Encoding

% The continuation problem encoded below includes seven monitor functions,
% and the corresponding inactive continuation parameters 'gap', 'p1', 'p2',
% 'eps1', 'eps2', 'y11e', and 'y21e'. Its dimensional deficit equals -1.
% The calls to the coco entry-point function indicate a desired manifold
% dimension of 1. To this end, two continuation parameters are released and
% allowed to vary during each continuation run.

p0   = [0.5; 0];
dev0 = [0.03; 0.2];
vu   = [sqrt(4*p0(1)+p0(2)^2)-p0(2); 2*p0(1)];
vu   = vu/norm(vu, 2);
vs   = [-sqrt(4*(1-p0(1))+p0(2)^2)-p0(2); 2*(1-p0(1))];
vs   = vs/norm(vs, 2);

% The initial trajectory segments may be discretized using a two-point time
% history, as in Recipes for Continuation.

% segs(1).t0 = [0; 1];
% x0         = dev0(1)*vu;
% segs(1).x0 = [x0  x0+huxley(x0, p0)]';
% segs(1).p0 = p0; 
% segs(2).t0 = [0; 1];
% x0         = [1; 0]+dev0(2)*vs;
% segs(2).x0 = [x0-huxley(x0, p0) x0]';
% segs(2).p0 = p0;

% Alternatively, we may use single-point time-histories

segs(1).t0 = 0;
segs(1).x0 = dev0(1)*vu';
segs(1).p0 = p0;
segs(2).t0 = 0;
segs(2).x0 = [1 0]+dev0(2)*vs';
segs(2).p0 = p0;

fprintf('\n Run=''%s'': Continue trajectory segments until y11e=0.5.\n', ...
  'huxley1');

prob = huxley_isol2het(coco_prob(), segs, dev0);
coco(prob, 'huxley1', [], 1, {'y11e', 'gap'}, [0 0.5]);

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s'' until y21e=0.5.\n', ...
  'huxley2', 5, 'huxley1');

prob = huxley_sol2het(coco_prob(), 'huxley1', 5);
coco(prob, 'huxley2', [], 1, {'y21e', 'gap'}, [0.5 1]);

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s'' until gap=0.\n', ...
  'huxley3', 2, 'huxley2');

prob = huxley_sol2het(coco_prob(), 'huxley2', 2);
coco(prob, 'huxley3', [], 1, {'gap', 'p2'}, [-0.2 0]);

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s'' until dev1=1e-3.\n', ...
  'huxley5', 4, 'huxley3');

prob = huxley_sol2het(coco_prob(), 'huxley3', 4);
coco(prob, 'huxley4', [], 1, {'dev1', 'p2'}, [1e-3 dev0(1)]);

fprintf(...
  '\n Run=''%s'': Continue segments from point %d in run ''%s'' until dev2=1e-3.\n', ...
  'huxley5', 3, 'huxley4');

prob = huxley_sol2het(coco_prob(), 'huxley4', 3);
coco(prob, 'huxley5', [], 1, {'dev2', 'p2'}, [1e-3 dev0(2)]);

fprintf(...
  '\n Run=''%s'': Continue constrained segments from point %d in run ''%s''.\n', ...
  'huxley6', 3, 'huxley5');

prob = huxley_sol2het(coco_prob(), 'huxley5', 3);
coco(prob, 'huxley6', [], 1, {'p1', 'p2'}, [0.25 0.75]);

% Use the following commands to graph a two-parameter bifurcation diagram
% corresponding to pairs of parameter values for points on the solution
% manifold. 

% bd6 = coco_bd_read('huxley6'); % Extract bifurcation data
% p1  = coco_bd_col(bd6, 'p1'); % Extract column data
% p2  = coco_bd_col(bd6, 'p2'); % Extract column data
% plot(p1, p2, 'r.', p1, (1-2*p1)/sqrt(2), 'k') % Compare with exact solution
%
% Alternatively, use the commands
% figure; hold on
% coco_plot_bd('huxley6','p1','p2')
% coco_plot_bd(struct('lspec',{{'ro', 'MarkerSize', 12}}), ...
%   'huxley6', 'p1', 'p1', @(p1) (1-2.*p1)/sqrt(2))
% hold off

%% Graphical representation of stored solutions

% Vector field
N     = 20;
x     = linspace(-0.1,1.1,2*N);
y     = linspace(-0.25,0.25,N);
[X, Y] = meshgrid(x,y);
fxy   = huxley([X(:) Y(:)]', repmat([0.5;0], [1 numel(X)]));
fx    = reshape(fxy(1,:), size(X));
fy    = reshape(fxy(2,:), size(X));

figure(1); clf; hold on; grid on; box on; axis([-0.1 1.1 -0.25 0.25])

quiver(X, Y, fx, fy, 1, 'LineStyle', '-', 'LineWidth',   1, ...
  'Color', [0.3 0.3 0.3], 'ShowArrowHead', 'on',  'MaxHeadSize', 0.2)
plot_skeleton(p0, vs, vu, 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black')
plot(0.5, 0, 'LineStyle', '-', 'LineWidth', sqrt(2), ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, ...
  'MarkerFaceColor', 'white')

xlabel('x_1'); ylabel('x_2')

hold off

%Stages of homotopy
labs = [1, 5, 2, 4, 3, 3];
runs = {'huxley1', 'huxley1', 'huxley2', 'huxley3', 'huxley4', 'huxley5'};

for i=1:6
  figure(i+1)
  clf
  hold on
  grid on
  box on
  axis([-0.1 1.1 -0.25 0.25])
  
  plot_skeleton(p0, vs, vu, 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', [0.7 0.7 0.7]);
  plot(0.5, 0, 'LineStyle', '-', 'LineWidth', sqrt(2), ...
    'Color', [0.5 0.5 0.5], 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'white');
  coco_plot_sol(runs{i}, labs(i), 'huxley', 1:2, 'x', 'x');
  
  hold off
end
