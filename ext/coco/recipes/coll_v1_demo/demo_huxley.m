coco_use_recipes_toolbox coll_v1 % add the coll_v1 toolbox to the search path

%% Section 7.3.2

% The continuation problem encoded below consists of a family of 202 zero
% functions (80 collocation conditions for each segment, 18 continuity
% conditions for each segment, two gluing conditions for duplication of
% problem parameters, and four boundary conditions) in terms of 208
% continuation variables (100 basepoint values for each segment, two
% interval lengths, two copies of the problem parameters, and two
% deviations of the heteroclinic trajectory end points from the
% corresponding equilibria), seven monitor functions, and the corresponding
% inactive continuation parameters 'gap', 'p1', 'p2', 'eps1', 'eps2',
% 'y11e', and 'y21e'. Its dimensional deficit equals -1. The calls to the
% coco entry-point function indicate a desired manifold dimension of 1. To
% this end, two continuation parameters are released and allowed to vary
% during each continuation run.

p0   = [0.5; 0];
eps0 = [0.03; 0.2];
vu   = [sqrt(4*p0(1)+p0(2)^2)-p0(2); 2*p0(1)];
vu   = vu/norm(vu, 2);
segs(1).t0 = [0; 1];
x0         = eps0(1)*vu;
segs(1).x0 = [x0  x0+huxley(x0, p0)]';
segs(1).p0 = p0;
vs   = [-sqrt(4*(1-p0(1))+p0(2)^2)-p0(2); 2*(1-p0(1))];
vs   = vs/norm(vs, 2);
segs(2).t0 = [0; 1];
x0         = [1; 0]+eps0(2)*vs;
segs(2).x0 = [x0-huxley(x0, p0) x0]';
segs(2).p0 = p0;

prob = huxley_isol2het(coco_prob(), segs, eps0);
coco(prob, 'huxley1', [], 1, {'y11e', 'gap'}, [0 0.5]);

prob = huxley_sol2het(coco_prob(), 'huxley1', 5);
coco(prob, 'huxley2', [], 1, {'y21e', 'gap'}, [0.5 1]);
prob = huxley_sol2het(coco_prob(), 'huxley2', 2);
coco(prob, 'huxley3', [], 1, {'gap', 'p2'}, [-0.2 0]);
prob = huxley_sol2het(coco_prob(), 'huxley3', 4);
coco(prob, 'huxley4', [], 1, {'eps1', 'p2'}, [1e-3 eps0(1)]);
prob = huxley_sol2het(coco_prob(), 'huxley4', 3);
coco(prob, 'huxley5', [], 1, {'eps2', 'p2'}, [1e-3 eps0(2)]);

% Use the following commands to graph the last discretized trajectory from
% run 'huxley5'. See also figure_7_3.m.

% hold on
% sol = coll_read_solution('huxley1', 'huxley5', 3); % Extract solution
% plot(sol.x(:,1), sol.x(:,2), 'r')
% sol = coll_read_solution('huxley2', 'huxley5', 3); % Extract solution
% plot(sol.x(:,1), sol.x(:,2), 'r')
% hold off

prob = huxley_sol2het(coco_prob(), 'huxley5', 3);
coco(prob, 'huxley6', [], 1, {'p1', 'p2'}, [0.25 0.75]);

% Use the following commands to graph a two-parameter bifurcation diagram
% corresponding to pairs of parameter values for points on the solution
% manifold. 

% bd6 = coco_bd_read('huxley6'); % Extract bifurcation data
% p1  = coco_bd_col(bd6, 'p1'); % Extract column data
% p2  = coco_bd_col(bd6, 'p2'); % Extract column data
% plot(p1, p2, 'r.', p1, (1-2*p1)/sqrt(2), 'k') % Compare with exact solution

coco_use_recipes_toolbox % remove the coll_v1 toolbox from the search path
