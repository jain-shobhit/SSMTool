coco_use_recipes_toolbox coll_v1 varcoll_v1 msbvp_v1 hspo_v2 atlas1d_v6 % Add atlas1d_v6 atlas algorithm and hspo_v2, msbvp_v1, varcoll_v1, and coll_v1 toolboxes to search path

%% Example 17.3

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
% respectively. By default, special points associated with saddle-node
% bifurcations, period-doubling bifurcations, and Neimark-Sacker
% bifurcations are detected by the 'hspo' toolbox and identified by the
% event types 'SN', 'PD', and 'NS'.

p0 = [0.2; 1; 1; 26; 1];
x0 = [0; 0; 0];
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

modes  = {'neg' 'pos'};
events = {'phase' 'phase'};
resets = {'phase' 'phase'};
t0     = {t1 t2};
x0     = {x1 x2};

prob = coco_set(coco_prob(), 'coll', 'NTST', 30);
prob = hspo_isol2segs(prob, '', ...
  {@duff, @duff_events, @duff_resets}, ...
  modes, events, resets, t0, x0, {'la' 'al' 'eps' 'A' 'om'}, p0);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'h', 2, 'PtMX', 200);
bd1 = coco(prob, 'run1', [], 1, 'om', [0.5 1.5]);

%% Example 17.5 [But variation in frequency, not amplitude]

% The continuation problem constructed below is identical to that above,
% but the initial solution guess is obtained from chart data stored in
% conjunction with the detection of a period-doubling bifurcation point in
% the previous run.

labs = coco_bd_labs(bd1, 'PD'); % Extract period-doubling solution labels
for lab=labs([1 end])
  [data chart] = coco_read_solution('', 'run1', lab); % Extract chart for first period-doubling point
  pd   = coco_get_chart_data(chart, 'hspo.PD'); % Extract initial solution guess for period-doubled branch
  
  prob = coco_set(coco_prob(), 'coll', 'NTST', 30);
  prob = hspo_isol2segs(prob, '', ...
    {@duff, @duff_events, @duff_resets}, ...
    pd.modes, pd.events, pd.resets, pd.t0, pd.x0, ...
    {'la' 'al' 'eps' 'A' 'om'}, p0);
  prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
  prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
  prob = coco_set(prob, 'cont', 'h', 2, 'PtMX', 200);
  bd2 = coco(prob, 'run2', [], 1, 'om', [0.5 1.5]);
  
  figure(1) 
  
  hold on
  om = coco_bd_col(bd2, 'om'); % Extract column data
  x0 = coco_bd_col(bd2, 'X0'); % Extract column data
  st = coco_bd_col(bd2,'hspo.test.stab'); % Extract column data
  plot(om(st==0),x0(1,st==0),'g.',om(st>0),x0(1,st>0),'r.')
  idx = coco_bd_idxs(bd2, 'EP'); % Extract row indices
  plot(om(idx),x0(1,idx),'go');
  idx = coco_bd_idxs(bd2, 'FP'); % Extract row indices
  plot(om(idx),x0(1,idx),'ko');
  idx = coco_bd_idxs(bd2, 'BP'); % Extract row indices
  plot(om(idx),x0(1,idx),'ro');
  idx = coco_bd_idxs(bd2, 'PD'); % Extract row indices
  plot(om(idx),x0(1,idx),'ro', 'MarkerFaceColor', 'y');
  hold off
  grid on
  drawnow
end

coco_use_recipes_toolbox % Remove atlas1d_v6 atlas algorithm and hspo_v2, msbvp_v1, varcoll_v1, and coll_v1 toolboxes from search path
