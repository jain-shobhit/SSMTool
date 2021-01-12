coco_use_recipes_toolbox hspo_v1 msbvp_v1 coll_v1 % add the hspo_v1, msbvp_v1, and coll_v1 toolboxes to the search path

%% Example 15.2

% The continuation problem structure encoded below corresponds to a family
% of 308 zero functions (120 collocation conditions per segment, 27
% continuity conditions per segment, eight boundary conditions, and six
% gluing conditions) in terms of 314 continuation variables (150 basepoint
% values per segment, two interval lengths, and two copies of the problem
% parameters), a family of six monitor functions that evaluate to the
% problem parameters and a single monitor function that evaluates to the
% second component (the physical velocity) of the end point at t=0 of the
% second trajectory segment, and the corresponding inactive continuation
% parameters 'k', 'c', 'A', 'w', 'd', and 'e', and active continuation
% parameter 'graze'. Its dimensional deficit equals 0. A one-dimensional
% family of impacting periodic orbits is obtained by releasing 'A' and
% allowing this to vary during continuation. Events associated with
% zero-crossings of the value of 'graze' (grazing periodic orbits) are
% detected and identified with the 'GR' event type.

% Continuation of impacting trajectory.
p0     = [1; 0.1; 1; 1; 1; 0.8];
modes  = {'free' 'free'};
events = {'impact' 'phase'};
resets = {'bounce' 'phase'};
f       = @(t, x) impact(x, p0, modes{1});
[t1 x1] = ode45(f, [0 3.2], [-0.98; -0.29; -pi]);
f       = @(t, x) impact(x, p0, modes{2});
[t2 x2] = ode45(f, [0 3.1], [1; -1.36; 0.076]);
t0 = {t1  t2};
x0 = {x1  x2};
hspo_args = {{@impact, @impact_events, @impact_resets}, ...
  modes, events, resets, t0, x0, {'k' 'c' 'A' 'w' 'd' 'e'}, p0};
prob = hspo_isol2segs(coco_prob(), '', hspo_args{:});

[data uidx] = coco_get_func_data(prob, 'msbvp.seg2.coll', ...
  'data', 'uidx'); % Extract toolbox data and context-dependent index array for the second segment
prob = coco_add_pars(prob, 'grazing', uidx(data.x0_idx(2)), ...
  'graze', 'active');
prob = coco_add_event(prob, 'GR', 'graze', 0);
prob = coco_set(prob, 'cont', 'ItMX', 100);
bd1 = coco(prob, 'impact1', [], {'A' 'graze'}, [0.01 1]);

% The continuation problem structure encoded below is identical to that
% above, but constructed from stored data associated with the grazing
% periodic orbit, and with a reassignment of 'graze' to the inactive
% continuation parameters and 'A' to the active continuation parameters. A
% one-dimensional family of grazing periodic orbits is obtained by
% releasing 'w' and allowing it to vary during continuation.

labgr = coco_bd_labs(bd1, 'GR');
prob = msbvp_sol2segs(coco_prob(), '', 'impact1', labgr);
[data uidx] = coco_get_func_data(prob, 'msbvp.seg2.coll', ...
  'data', 'uidx');
prob = coco_add_pars(prob, 'grazing', uidx(data.x0_idx(2)), ...
  'graze', 'active');
prob = coco_xchg_pars(prob, 'graze', 'A');
prob = coco_set(prob, 'cont', 'ItMX', 100);
coco(prob, 'impact2', [], {'w' 'A' 'graze'}, {[] [0 1]});

coco_use_recipes_toolbox % remove the hspo_v1, msbvp_v1, and coll_v1 toolboxes from the search path
