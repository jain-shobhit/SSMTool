coco_use_recipes_toolbox hspo_v1 msbvp_v1 coll_v1 % add the hspo_v1, msbvp_v1, and coll_v1 toolboxes to the search path

%% Example 9.3

% The continuation problem structure encoded below corresponds to a family
% of 307 zero functions (120 collocation conditions per segment, 27
% continuity conditions per segment, eight boundary conditions, and five
% gluing conditions) in terms of 312 continuation variables (150 basepoint
% values per segment, two interval lengths, and two copies of the problem
% parameters), a family of five monitor functions that evaluate to the
% problem parameters, and the corresponding inactive continuation
% parameters 'V', 'c', 'n', 'w', and 'e'. Its dimensional deficit equals 0.
% A one-dimensional family of non-impacting sticking periodic orbits is
% obtained by releasing 'V' and allowing this to vary during continuation.

p0 = [0.59; 0.04; 3.17; 0.92; 0.80];
modes  = {'stick' 'stick'};
events = {'phase' 'minsep'};
resets = {'phase' 'turn'};
f  = @(t, x) stickslip(x, p0, modes{1});
[t1, x1] = ode45(f, [0  1.5], [0.5; 0; 0]);
f  = @(t, x) stickslip(x, p0, modes{2});
[t2, x2] = ode45(f, [0  1.5], [0.25; 0; -pi/2]);
t0 = {t1 t2};
x0 = {x1 x2};
prob = hspo_isol2segs(coco_prob(), '', ...
  {@stickslip, @stickslip_events, @stickslip_resets}, ...
  {@stickslip_DFDX, [], []}, modes, events, resets, ...
  t0, x0, {'V' 'c' 'n' 'w' 'e'}, p0);
coco(prob, 'stickslip1', [], 1, 'V', [0.5 0.7]);

% The continuation problem structure encoded below includes an additional
% monitor function that evaluates to the first component of the end point
% at t=0 of the first trajectory segment, with corresponding inactive
% continuation parameter 'pos', whose initial value is assigned to equal
% 0.5. Its dimensional deficit equals -1. A unique solution to a closed
% continuation problem is obtained by releasing 'V', while keeping 'pos'
% fixed. The corresponding periodic orbit grazes the impact surface.

[data uidx] = coco_get_func_data(prob, 'msbvp', 'data', 'uidx');  % Extract 'msbvp' toolbox data and context-dependent index set
prob = coco_add_pars(prob, 'graze', uidx(data.x0_idx(1)), 'pos');
prob = coco_set_parival(prob, 'pos', 0.5);
coco(prob, 'graze_run', [], 0, {'V' 'pos'});

coco_use_recipes_toolbox % remove the hspo_v1, msbvp_v1, and coll_v1 toolboxes from the search path

coco_use_recipes_toolbox hspo_v1 msbvp_v1 coll_v2 % add the hspo_v1, msbvp_v1, and coll_v2 toolboxes to the search path

%% Example 9.4

% The continuation problem structure encoded below corresponds to a family
% of 513 zero functions (120 collocation conditions and 27 continuity
% conditions each for the first and second segments, 160 collocation
% conditions and 36 continuity conditions for the third segment, 13
% boundary conditions, and 10 gluing conditions) in terms of 518
% continuation variables (150 basepoint values each for the first and
% second segments, 200 basepoint values for the third segment, three
% interval lengths, and three copies of the problem parameters), a family
% of five monitor functions that evaluate to the problem parameters, and
% the corresponding inactive continuation parameters 'V', 'c', 'n', 'w',
% and 'e'. A one-dimensional family of periodic, impacting stick-slip
% oscillations results by releasing 'V' and allowing it to vary during
% continuation.

sol = msbvp_read_solution('', 'graze_run', 1);
modes  = {'stick' 'stick' 'slip'};
events = {'phase' 'collision' 'rest'};
resets = {'phase' 'bounce' 'stick'};
t0 = {sol{1}.t, sol{2}.t, 0};                 % Third segment has zero length
x0 = {sol{1}.x, sol{2}.x, [0 sol{1}.x(1,:)]}; % And consists of a single point
p0 = sol{1}.p; % This line is missing from Recipes for Continuation
prob = hspo_isol2segs(coco_prob(), '', ...
  {@stickslip, @stickslip_events, @stickslip_resets}, ...
  {@stickslip_DFDX, [], []}, modes, events, resets, ...
  t0, x0, {'V' 'c' 'n' 'w' 'e'}, p0);
coco(prob, 'stickslip2', [], 1, 'V', [0.5 0.7]);

coco_use_recipes_toolbox % remove the hspo_v1, msbvp_v1, and coll_v2 toolboxes from the search path
