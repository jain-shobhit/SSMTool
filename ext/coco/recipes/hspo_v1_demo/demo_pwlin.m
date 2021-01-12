coco_use_recipes_toolbox hspo_v1 msbvp_v1 coll_v1 % add the hspo_v1, msbvp_v1, and coll_v1 toolboxes to the search path

%% Example 9.2

% The continuation problem structure encoded below corresponds to a family
% of 345 zero functions (120 collocation conditions for the first segment,
% 18 continuity conditions for the first segment, 160 collocation
% conditions for the second segment, 38 continuity conditions for the
% second segment, six boundary conditions, and three gluing conditions) in
% terms of 348 continuation variables (140 basepoint values for the first
% segment, 200 basepoint values for the second segment, two interval
% lengths, and two copies of the problem parameters), three monitor
% functions that evaluate to the problem parameters, and the corresponding
% inactive continuation parameters 'al', 'be', and 'ga'. Its dimensional
% deficit equals 0. To obtain a one-dimensional family of periodic orbits,
% the continuation parameter 'be' is released and allowed to vary during
% continuation.

p0     = [1; 2; 1.5];
modes  = {'left' 'right'};
events = {'boundary' 'boundary'};
resets = {'switch' 'switch'};
t0 = linspace(0, pi, 100)';
x1 = [-sin(t0) 0.5+cos(t0)];
x2 = [ sin(t0) 0.5-cos(t0)];
t0 = {t0 0.5*t0};
x0 = {x1 x2};

prob = coco_prob();
prob = coco_set(prob, 'msbvp.seg1.coll', 'NTST', 10, 'NCOL', 6);
prob = coco_set(prob, 'msbvp.seg2.coll', 'NTST', 20, 'NCOL', 4);
prob = hspo_isol2segs(prob, '', ...
  {@pwlin, @pwlin_events, @pwlin_resets}, ...
  {@pwlin_DFDX, @pwlin_events_DFDX, @pwlin_resets_DFDX}, ...
  {@pwlin_DFDP, @pwlin_events_DFDP, @pwlin_resets_DFDP}, ...
  modes, events, resets, t0, x0, {'al' 'be' 'ga'}, p0);
coco(prob, 'pwlin1', [], 1, 'be', [0 5]);

% The continuation problem encoded below is identical to that above, but
% constructed from stored data obtained in the previous run. Toolbox
% settings encoded in the stored toolbox data override the default values.
% A one-dimensional family of periodic orbits is obtained by releasing 'al'
% and allowing it to vary during continuation.

prob = msbvp_sol2segs(coco_prob(), '', 'pwlin1', 9);
coco(prob, 'pwlin2', [], 1, 'al', [0 4]);

coco_use_recipes_toolbox % remove the hspo_v1, msbvp_v1, and coll_v1 toolboxes from the search path
