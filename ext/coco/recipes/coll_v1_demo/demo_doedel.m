coco_use_recipes_toolbox coll_v1 alg_v6 % add the coll_v1 and alg_v6 toolboxes to the search path

%% Section 7.3.3

% The continuation problem encoded below consists of a family of 213 zero
% functions (80 collocation conditions for each segment, 18 continuity
% conditions for each segment, two algebraic equations for each
% equilibrium, six gluing conditions for duplication of problem parameters,
% three eigenspace conditions, and four boundary conditions) in terms of
% 220 continuation variables (100 basepoint values for each segment, two
% interval lengths, two problem variables for each equilibrium, four copies
% of the problem parameters, two components of an eigenvector and the
% corresponding eigenvalue, two deviations of the heteroclinic trajectory
% end points from the corresponding equilibria, and the angle with the
% horizontal), seven monitor functions, and the corresponding inactive
% continuation parameters 'gap', 'p1', 'p2', 'eps1', 'eps2', 'y12e', and
% 'y22e'. Its dimensional deficit equals 0. The calls to the coco
% entry-point function indicate a desired manifold dimension of 1. To this
% end, one continuation parameters is released and allowed to vary during
% each continuation run.

p0 = [1; 1];
eps0 = [0.05; 0.05];
th0 = -pi/2;
eqs10 = [-1; 1]; % Equilibrium
eqs20 = [1; -1]; % Equilibrium
vec0 = [-3/sqrt(10); 1/sqrt(10)];
lam0 = -2;
segs(1).t0 = [0; 1];
x0         = eqs10+eps0(1)*[cos(th0); sin(th0)];
segs(1).x0 = [x0  x0+doedel(x0, p0)]'; % Euler step
segs(1).p0 = p0;
segs(2).t0 = [0; 1];
x0         = eqs20+eps0(2)*vec0;
segs(2).x0 = [x0-doedel(x0, p0) x0]'; % Euler step
segs(2).p0 = p0;
algs(1).x0 = eqs10;
algs(1).p0 = p0;
algs(2).x0 = eqs20;
algs(2).p0 = p0;

prob = coco_prob();
prob = doedel_isol2het(prob, segs, algs, eps0, th0, vec0, lam0);
coco(prob, 'doedel1', [], 1, 'y12e', [0 0.99]);
prob = doedel_sol2het(coco_prob(), 'doedel1', 3);
coco(prob, 'doedel2', [], 1, 'y22e', [-0.995 0]);
prob = doedel_sol2het(coco_prob(), 'doedel2', 5);
coco(prob, 'doedel3', [], 1, 'gap', [-2 0]);
prob = doedel_sol2het(coco_prob(), 'doedel3', 5);
coco(prob, 'doedel4', [], 1, 'eps1', [1e-3 eps0(1)]);
prob = doedel_sol2het(coco_prob(), 'doedel4', 3);
coco(prob, 'doedel5', [], 1, 'eps2', [1e-3 eps0(2)]);
prob = doedel_sol2het(coco_prob(), 'doedel5', 2);
coco(prob, 'doedel6', [], 1, 'p2', [0.5 8]);

coco_use_recipes_toolbox % remove the coll_v1 and alg_v6 toolboxes from the search path
