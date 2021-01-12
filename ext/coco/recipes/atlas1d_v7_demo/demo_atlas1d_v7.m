%% Example 17.4 [No code shown in Recipes for Continuation]

% The continuation problem structure encoded below corresponds to a family
% of 907 zero functions (360 collocation conditions per segment, 87
% continuity conditions per segment, 8 boundary conditions, and five gluing
% conditions) in terms of 912 continuation variables (450 basepoint values
% per segment, two interval lengths, and two pairs of parameter values), a
% family of five monitor functions that evaluate to the problem parameters,
% and the corresponding inactive continuation parameters 'la', 'al', 'eps',
% 'A', and 'om'. Its dimensional deficit equals 0. A one-dimensional family
% of periodic orbits is obtained by releasing 'A' and allowing it vary
% during continuation. The atlas algorithm detects fold points and branch
% points and identifies these by the event types 'FP' and 'BP',
% respectively. 

%% Without automatic branch switching.

coco_use_recipes_toolbox coll_v1 msbvp_v1 hspo_v2 varcoll_v1 atlas1d_v6 % Add atlas1d_v6 atlas algorithm and hspo_v1, msbvp_v1, and coll_v1 toolboxes to search path

p0 = [0.2 1 1 1 1]';
x0 = [0 0 0]';
for i=1:10 % Let transients die out
  [t0 x0] = ode45(@(t,x) duff(x,p0,'neg'), [0 pi], x0); %#ok<*ASGLU>
  x0 = x0(end,:)';
  [t0 x0] = ode45(@(t,x) duff(x,p0,'pos'), [0 pi], x0);
  x0 = x0(end,:)';
end

x0 = [x0(1:2);0];
[t1 x1] = ode45(@(t,x) duff(x,p0,'neg'), [0 pi], x0); % First segment of approximate periodic orbit
x0 =  [x1(end,1:2) 0]';
[t2 x2] = ode45(@(t,x) duff(x,p0,'pos'), [0 pi], x0); % Second segment of approximate periodic orbit

modes  = {'neg'   'pos'  };
events = {'phase' 'phase'};
resets = {'phase' 'phase'};
t0     = {  t1,     t2};
x0     = {  x1,     x2};

prob = coco_set(coco_prob(), 'coll', 'NTST', 30);
prob = hspo_isol2segs(prob, '', ...
  {@duff, @duff_events, @duff_resets}, ...
  modes, events, resets, t0, x0, ...
  {'la' 'al' 'eps' 'A' 'om'}, p0);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat'); % Store initial point to bifurcation data
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'h', 1, 'almax', 30, 'PtMX', 5000);
coco(prob, 'run1', [], 1, 'A', [0.5 81]);

coco_use_recipes_toolbox % Remove atlas1d_v6 atlas algorithm and hspo_v2, msbvp_v1, and coll_v1 toolboxes from search path

%% With automatic branch switching.

coco_use_recipes_toolbox coll_v1 msbvp_v1 hspo_v2 varcoll_v1 atlas1d_v7 % Add atlas1d_v7 atlas algorithm and hspo_v1, msbvp_v1, and coll_v1 toolboxes to search path

p0 = [0.2 1 1 1 1]';
x0 = [0 0 0]';
for i=1:10 % Let transients die out
  [t0 x0] = ode45(@(t,x) duff(x,p0,'neg'), [0 pi], x0);
  x0 = x0(end,:)';
  [t0 x0] = ode45(@(t,x) duff(x,p0,'pos'), [0 pi], x0);
  x0 = x0(end,:)';
end

x0 = [x0(1:2);0];
[t1 x1] = ode45(@(t,x) duff(x,p0,'neg'), [0 pi], x0); % First segment of approximate periodic orbit
x0 =  [x1(end,1:2) 0]';
[t2 x2] = ode45(@(t,x) duff(x,p0,'pos'), [0 pi], x0); % Second segment of approximate periodic orbit

modes  = {'neg'   'pos'  };
events = {'phase' 'phase'};
resets = {'phase' 'phase'};
t0     = {  t1,     t2   };
x0     = {  x1,     x2   };

prob = coco_set(coco_prob(), 'coll', 'NTST', 30);
prob = hspo_isol2segs(prob, '', ...
  {@duff, @duff_events, @duff_resets}, ...
  modes, events, resets, t0, x0, ...
  {'la' 'al' 'eps' 'A' 'om'}, p0);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat'); % Store initial point to bifurcation data
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'h', 1, 'almax', 30, 'PtMX', 5000);
coco(prob, 'run2', [], 1, 'A', [0.5 81]);

coco_use_recipes_toolbox % Remove atlas1d_v6 atlas algorithm and hspo_v1, msbvp_v1, and coll_v1 toolboxes from search path

%% Compute branch of period-2 orbits

coco_use_recipes_toolbox coll_v1 msbvp_v1 hspo_v2 varcoll_v1 atlas1d_v6 % Add atlas1d_v6 atlas algorithm and hspo_v1, msbvp_v1, and coll_v1 toolboxes to search path

rrunid = 'run2';
bd = coco_bd_read(rrunid);
labs = coco_bd_labs(bd, 'PD');
rlab = labs(end);

[data chart] = coco_read_solution('', rrunid, rlab);
pd   = coco_get_chart_data(chart, 'hspo.PD');

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30);
prob = hspo_isol2segs(prob, '', ...
  {@duff, @duff_events, @duff_resets}, ...
  pd.modes, pd.events, pd.resets, pd.t0, pd.x0, ...
  {'la' 'al' 'eps' 'A' 'om'}, p0);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'h', 1, 'almax', 30, 'PtMX', 5000);
coco(prob, 'run3', [], 1, 'A', [0.5 81]);

%% Compute branch of period-4 orbits

rrunid = 'run3';
bd = coco_bd_read(rrunid);
labs = coco_bd_labs(bd, 'PD');
rlab = labs(end);

[data chart] = coco_read_solution('', rrunid, rlab);
pd   = coco_get_chart_data(chart, 'hspo.PD');

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30);
prob = hspo_isol2segs(prob, '', ...
  {@duff, @duff_events, @duff_resets}, ...
  pd.modes, pd.events, pd.resets, pd.t0, pd.x0, ...
  {'la' 'al' 'eps' 'A' 'om'}, p0);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'h', 2, 'almax', 30, 'PtMX', 5000);
coco(prob, 'run4', [], 1, 'A', [0.5 81]);

%% Compute branch of period-8 orbits

rrunid = 'run4';
bd = coco_bd_read(rrunid);
labs = coco_bd_labs(bd, 'PD');
rlab = labs(end);

[data chart] = coco_read_solution('', rrunid, rlab);
pd   = coco_get_chart_data(chart, 'hspo.PD');

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30);
prob = hspo_isol2segs(prob, '', ...
  {@duff, @duff_events, @duff_resets}, ...
  pd.modes, pd.events, pd.resets, pd.t0, pd.x0, ...
  {'la' 'al' 'eps' 'A' 'om'}, p0);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'h', 4, 'almax', 30, 'PtMX', 5000);
bd1 = coco(prob, 'run5', [], 1, 'A', [0.5 81]);

coco_use_recipes_toolbox % Remove atlas1d_v6 atlas algorithm and hspo_v1, msbvp_v1, and coll_v1 toolboxes from search path
