coco_use_recipes_toolbox po_v3 coll_v6 % add the coll_v6 and po_v3 toolboxes to the search path

%% Section 20.2.2

% The continuation problem structure encoded below initially consists of
% 101 zero functions (80 collocation conditions, 18 continuity conditions,
% two periodicity conditions, and one integral phase condition) in terms of
% 102 continuation variables (100 basepoint values, one interval length,
% and one problem parameter), a family of two monitor functions that
% evaluate to the interval length and the problem parameter, respectively,
% with the corresponding active continuation parameter 'po.period' and
% inactive continuation parameter 'eps'. Its dimensional deficit equals 0.
% A family of approximate periodic orbits is obtained by releasing 'eps'
% and allowing it to vary during continuation. Terminal special points
% associated with a discretization error estimate exceeding a given
% tolerance (here 1e-3) are identified by 'MCXL'. Adaptive remeshing is
% applied to the continuation problem after every step, corresponding to an
% almost equidistributed the temporal mesh with respect to an error
% estimate.

eps0 = 0.1;
t0 = linspace(0,2*pi,100)';
x0 = [ sin(t0) cos(t0) ];

prob  = coco_prob();
prob  = coco_set(prob, 'cont', 'ItMX', 100);
prob  = coco_set(prob, 'cont', 'NAdapt', 1); % Remesh after every continuation step
prob  = coco_set(prob, 'coll', 'TOL', 1.0e-3);

prob2 = po_isol2orb(prob, '', @pneta, t0, x0, 'eps', eps0);
bd = coco(prob2, 'run1', [], 1, {'eps' 'po.period' 'po.seg.coll.err' 'po.seg.coll.NTST'}, [0.1 20]);

% The continuation problem encoded below is constructed using the last
% solution obtained in the previous run with discretization error estimate
% below the given tolerance. It initialy consists of 311 zero functions
% (248 collocation conditions, 60 continuity conditions, two periodicity
% conditions, and one integral phase condition) in terms of 312
% continuation variables (310 basepoint values, one interval length, and
% one problem parameter), a family of two monitor functions that evaluate
% to the interval length and the problem parameter, respectively, with the
% corresponding active continuation parameter 'po.period' and inactive
% continuation parameter 'eps'. Its dimensional deficit again equals 0. A
% family of approximate periodic orbits is obtained by releasing 'eps' and
% allowing it to vary during continuation. Terminal special points
% associated with a discretization error estimate exceeding a given
% tolerance (here 1e-3) are identified by 'MCXL'.Adaptive remeshing is
% applied to the continuation problem after every step, corresponding to an
% almost equidistributed the temporal mesh with respect to an error
% estimate.

labs = coco_bd_labs(bd, 'EP');
prob2 = po_sol2orb(prob, '', 'run1', labs(end));
coco(prob2, 'run2', [], 1, {'eps' 'po.seg.coll.err' 'po.seg.coll.NTST'}, [0.1 20]);

coco_use_recipes_toolbox % Remove the coll_v6 and po_v3 toolboxes from the search path
