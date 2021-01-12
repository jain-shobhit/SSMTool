coco_use_recipes_toolbox po_v3 coll_v4 % add the coll_v4 and po_v3 toolboxes to the search path

%% A comoving mesh [Code not in Recipes for Continuation]

% The continuation problem structure encoded below consists of 101 zero
% functions (80 collocation conditions, 18 continuity conditions, two
% periodicity conditions, and one integral phase condition) in terms of 102
% continuation variables (100 basepoint values, one interval length, and
% one problem parameter), a family of two monitor functions that evaluate
% to the interval length and the problem parameter, respectively, with the
% corresponding active continuation parameter 'po.period' and inactive
% continuation parameter 'eps'. Its dimensional deficit equals 0. A family
% of approximate periodic orbits is obtained by releasing 'eps' and
% allowing it to vary during continuation. Terminal special points
% associated with a discretization error estimate exceeding a given
% tolerance (here 1e-3) are identified by 'MCXL'.

eps0 = 0.1;
t0   = linspace(0, 2*pi, 100)';
x0   = [sin(t0) cos(t0)];
prob  = coco_set(coco_prob(), 'cont', 'ItMX', [0 200]);
prob  = coco_set(prob, 'coll', 'TOL', 1.0e-3);

prob2 = po_isol2orb(prob, '', @pneta, t0, x0, 'eps', eps0);
coco(prob2, 'run1', [], 1, {'eps' 'po.seg.coll.err'}, [0.1 20]);

% The continuation problem encoded below is constructed using the last
% solution obtained in the previous run with discretization error estimate
% below the given tolerance. It consists of 201 zero functions (160
% collocation conditions, 38 continuity conditions, two periodicity
% conditions, and one integral phase condition) in terms of 202
% continuation variables (200 basepoint values, one interval length, and
% one problem parameter), a family of two monitor functions that evaluate
% to the interval length and the problem parameter, respectively, with the
% corresponding active continuation parameter 'po.period' and inactive
% continuation parameter 'eps'. Its dimensional deficit again equals 0. A
% family of approximate periodic orbits is obtained by releasing 'eps' and
% allowing it to vary during continuation. Terminal special points
% associated with a discretization error estimate exceeding a given
% tolerance (here 1e-3) are identified by 'MCXL'.

bd  = coco_bd_read('run1'); % Extract bifurcation data
lab = coco_bd_labs(bd, 'MXCL'); % Extract MXCL solution label
prob  = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob2 = po_sol2orb(prob, '', 'run1', lab);
coco(prob2, 'run2', [], 1, {'eps' 'po.seg.coll.err'}, [0.1 20]);

% The continuation problem encoded below is constructed using the last
% solution obtained in the previous run with discretization error estimate
% below the given tolerance. It consists of 601 zero functions (500
% collocation conditions, 98 continuity conditions, two periodicity
% conditions, and one integral phase condition) in terms of 602
% continuation variables (600 basepoint values, one interval length, and
% one problem parameter), a family of two monitor functions that evaluate
% to the interval length and the problem parameter, respectively, with the
% corresponding active continuation parameter 'po.period' and inactive
% continuation parameter 'eps'. Its dimensional deficit again equals 0. A
% family of approximate periodic orbits is obtained by releasing 'eps' and
% allowing it to vary during continuation. Terminal special points
% associated with a discretization error estimate exceeding a given
% tolerance (here 1e-3) are identified by 'MCXL'.

bd  = coco_bd_read('run2'); % Extract bifurcation data
lab = coco_bd_labs(bd, 'MXCL'); % Extract MXCL solution label
prob  = coco_set(prob, 'coll', 'NTST', 50, 'NCOL', 5);
prob2 = po_sol2orb(prob, '', 'run2', lab);
coco(prob2, 'run3', [], 1, {'eps' 'po.seg.coll.err'}, [0.1 20]);

coco_use_recipes_toolbox % remove the coll_v4 and po_v3 toolboxes from the search path
