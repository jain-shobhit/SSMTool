coco_use_recipes_toolbox varcoll_v2 po_v1 coll_v2 % Add po_v1, varcoll_v2, and coll_v2 toolboxes to search path

%% Stage I

% The continuation problem constructed below corresponds to a family of 601
% zero functions (120 collocation conditions, 27 continuity conditions,
% three periodicity conditions, one integral phase condition, 360
% variational collocation conditions, 81 variational continuity conditions,
% and 9 boundary conditions) in terms of 604 continuation variables (150
% trajectory basepoint values, one interval length, three problem
% parameters, and 450 fundamental solution basepoint values), a family of
% four monitor functions that evaluate to the problem parameters and the
% period, and the corresponding inactive continuation parameters, 's', 'b',
% and 'r', and active continuation parameter 'po.period'. Its dimensional
% deficit equals 0. A family of periodic orbits emanating from a Hopf
% bifurcation is obtained by releasing 'r' and allowing it to vary during
% continuation.

s  = 10;
b  = 8/3;
r  = 470/19;
p0 = [s; r; b];
eq = [-sqrt(b*(r-1)) -sqrt(b*(r-1)) r-1];       % Equilibrium branch
om = 4*sqrt(110/19);                            % Angular frequency at Hopf bifurcation
re = [-20/9*sqrt(38/1353) 2/9*sqrt(38/1353) 1]; % Real part of eigenvector
im = [-19/9*sqrt(5/123) -35/9*sqrt(5/123) 0];   % Imaginary part of eigenvector
t0 = (0:2*pi/100:2*pi)'/om;                                   % Temporal mesh
x0 = repmat(eq, size(t0))+0.01*(cos(om*t0)*re-sin(om*t0)*im); % Trajectory on mesh
coll_args = {@lorenz, @lorenz_DFDX, @lorenz_DFDP, t0, x0, ...
  {'s' 'r' 'b'}, p0}; % [Typo in Recipes for Continuation, 1st edition: Lorenz, not Lorentz!]
var_args = {@lorenz_DFDXDX, @lorenz_DFDXDP};
prob = povar_isol2orb(coco_prob(), '', coll_args{:}, var_args{:}); % Construct orbit problem with variational equations
coco(prob, 'runHopf', [], 1, 'r', [24 25]);

%% Stage II

% The continuation problem structure below consists of an additional 310
% zero functions (120 collocation conditions each for two segments, 27
% continuity conditions each for two segments, two gluing conditions for
% redundant copies of the problem parameters, four eigenvector conditions,
% and six eigenspace conditions) in terms of an additional 314 continuation
% variables (150 basepoint values each for two segments, two interval
% lengths, two additional copies of the problem parameters, an eigenvector,
% a Floquet multiplier, and two separations from end points to orbit and
% equilibrium), an additional family of six monitor functions that evaluate
% to two hyperplane projections and the end point separations and interval
% lenghts, respectively, and the corresponding inactive continuation
% parameters 'sg1', 'sg2', 'eps1', 'eps2', 'T1', 'T2'. Its dimensional
% deficit equals -2. We grow an orbit in the unstable manifold of the
% equilibrium by releasing 'sg1', 'T1', and 'T2', and allow these to vary
% during continuation.

prob = riess_start_1(coco_prob(), 'runHopf', 6);
prob = coco_set(prob, 'cont', 'ItMX', 500);
prob = coco_set(prob, 'cont', 'NPR', 50);
cont_args = {{'sg1' 'T1'  'T2'}, {[-30 0] [0 1]}};
coco(prob, 'run1', [], 1, cont_args{:});

%% Stage III

% The continuation problem structure encoded below is identical to that
% above, but constructed from stored data. We grow an orbit in the stable
% manifold of the periodic orbit by releasing 'sg2', 'T2', and 'T1', and
% allow these to vary during continuation.

prob = riess_restart_1(coco_prob(), 'run1', 6);
prob = coco_set(prob, 'cont', 'ItMX', 500);
prob = coco_set(prob, 'cont', 'NPR', 50);
cont_args = {{'sg2' 'T2' 'T1'}, {[0 30] [0 1]}};
coco(prob, 'run2', [], 1, cont_args{:});

%% Stage IV

% The continuation problem structure encoded below is identical to that
% above, but constructed from stored data. We sweep a family of orbits in the stable
% manifold of the periodic orbit by releasing 'eps2', 'T1', and 'T2', and
% allow these to vary during continuation.

prob = riess_restart_1(coco_prob(), 'run2', 9);
prob = coco_set(prob, 'cont', 'ItMX', 300);
cont_args = {{'eps2' 'T1' 'T2'}, [1e-6 1e-1]};
coco(prob, 'run3', [], 1, cont_args{:});

%% Stage V

% The continuation problem structure below consists of an additional one
% zero function, an additional monitor function, and the corresponding
% inactive continuation parameter 'lingap'. Its dimensional deficit equals
% -4. We reduce the Lin gap to 0 by releasing 'lingap', 'r', 'eps2', 'T1',
% and 'T2', and allow these to vary during continuation.

prob = riess_start_2(coco_prob(), 'run3');
prob = coco_set(prob, 'cont', 'ItMX', 500);
cont_args = {{'lingap' 'r' 'eps2' 'T1' 'T2'}, [-1 0]};
coco(prob, 'run4', [], 1, cont_args{:});

%% Stage VI

% The continuation problem structure encoded below is identical to that
% above, but constructed from stored data. We continue in the problem
% parameters by releasing 'r', 'b', 'eps1', 'eps2', and 'T2', and allow
% these to vary during continuation.

% Continue in problem parameters
prob = riess_restart_2(coco_prob(), 'run4', 7);
prob = coco_set(prob, 'cont', 'ItMX', 500);
prob = coco_set(prob, 'cont', 'NPR', 50);
cont_args = {{'r' 'b' 'eps1' 'eps2' 'T2'}, [20 30]};
coco(prob, 'run5', [], 1, cont_args{:});

coco_use_recipes_toolbox % Remove po_v1, varcoll_v2, and coll_v2 toolboxes from search path
