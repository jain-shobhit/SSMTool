coco_use_recipes_toolbox po_v1 coll_v1 % add the coll_v1 and po_v1 toolboxes to the search path

%% Example 18.1

% The continuation problem structure encoded below consists of 101 zero
% functions (80 collocation conditions, 18 continuity conditions, two
% periodicity conditions, and one integral phase condition) in terms of 102
% continuation variables (100 basepoint values, one interval length, and
% one problem parameter), a family of two monitor functions that evaluate
% to the interval length and the problem parameter, respectively, with the
% corresponding active continuation parameter 'po.period' and inactive
% continuation parameter 'eps'. Its dimensional deficit equals 0. A family
% of approximate periodic orbits is obtained by releasing 'eps' and
% allowing it to vary during continuation.

eps0  = 0.1;
t0    = linspace(0, 2*pi, 100)';
x0    = [sin(t0) cos(t0)];
prob  = coco_set(coco_prob(), 'cont', 'ItMX', [0 200]);
prob1 = po_isol2orb(prob, '', @pneta, t0, x0, 'eps', eps0);
coco(prob1, 'pneta1', [], 1, 'eps', [0.1 20]);

% The continuation problem structure encoded below consists of 1801 zero
% functions (1500 collocation conditions, 298 continuity conditions, two
% periodicity conditions, and one integral phase condition) in terms of
% 1802 continuation variables (1800 basepoint values, one interval length,
% and one problem parameter), a family of two monitor functions that
% evaluate to the interval length and the problem parameter, respectively,
% with the corresponding active continuation parameter 'po.period' and
% inactive continuation parameter 'eps'. Its dimensional deficit equals 0.
% A family of approximate periodic orbits is againobtained by releasing
% 'eps' and allowing it to vary during continuation.

prob2 = coco_set(prob, 'coll', 'NTST', 150, 'NCOL', 5);
prob2 = po_isol2orb(prob2, '', @pneta, t0, x0, 'eps', eps0);
coco(prob2, 'pneta2', [], 1, 'eps', [0.1 20]);

coco_use_recipes_toolbox % remove the coll_v1 and po_v1 toolboxes from the search path
