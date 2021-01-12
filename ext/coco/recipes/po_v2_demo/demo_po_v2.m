coco_use_recipes_toolbox po_v2 varcoll_v1 coll_v1 % add the coll_v1, varcoll_v1, and po_v2 toolboxes to the search path

%% Detecting bifurcations and branch switching at period doubling

% The continuation problem constructed below consists of 151 zero functions
% (120 collocation conditions, 27 continuity conditions, three periodicity
% conditions, and one integral phase condition) in terms of 157
% continuation variables (150 basepoint values, one interval length, and
% six problem parameters), a family of seven monitor functions that
% evaluate to the problem parameters and the interval length, respectively,
% and the corresponding inactive continuation parameters 'nu', 'be', 'ga',
% 'r', 'a3', and 'b3' and active continuation parameter 'po.period'. Its
% dimensional deficit equals 0. A one-dimensional manifold of periodic
% orbits is obtained by releasing 'nu' and allowing it to vary during
% continuation. Special points associated with saddle-node bifurcations,
% period-doubling bifurcations, and Neimark-Sacker bifurcations are
% detected during continuation and denoted by the event types 'SN', 'PD',
% and 'NS'.

p0 = [-0.59 ; 0.5 ; -0.6 ; 0.6 ; 0.328578 ; 0.933578];
x0 = [0.177 ; 0.082 ; 0.183];
T0 = 9;
N  = 1;
[t0 x0] = ode45(@(t,x) tor(x,p0), [0 T0/N], x0(:,end)); % Approximate periodic orbit

prob = po_isol2orb(coco_prob(), '', @tor, @tor_DFDX, @tor_DFDP, t0, x0, ...
  {'nu' 'be' 'ga' 'r' 'a3' 'b3'}, p0);
prob = coco_set(prob, 'cont', 'FP', false, 'BP', false); % Turn off fold and branch point detection by atlas
bd1 = coco(prob, 'run1', [], 1, {'nu' 'po.period'}, [-0.65, -0.55]);

% The continuation problem constructed below is identical to that above,
% but the initial solution guess is obtained from chart data stored in
% conjunction with the detection of a period-doubling bifurcation point in
% the previous run.

labs = coco_bd_labs(bd1, 'PD'); % Extract labels for period-doubling bifurcation
[data chart] = coco_read_solution('', 'run1', labs(1)); % Extract chart for first period-doubling point
pd   = coco_get_chart_data(chart, 'po.PD'); % Extract initial solution guess for period-doubled branch
prob = po_isol2orb(coco_prob(), '', @tor, @tor_DFDX, @tor_DFDP, pd.pd_t0, ...
  pd.pd_x0, {'nu' 'be' 'ga' 'r' 'a3' 'b3'}, pd.pd_p);
prob = coco_set(prob, 'cont', 'FP', false, 'BP', false);
bd2 = coco(prob, 'run2', [], 1, {'nu' 'po.period'}, [-1, -0.55]);

labs = coco_bd_labs(bd2, 'all'); % Extract solution labels
clf
hold on
grid on
for lab=labs
  sol = po_read_solution('','run2', lab); % Extract trajectory segment
  plot3(sol.x(:,1),sol.x(:,2),sol.x(:,3),'r')
end
hold off

coco_use_recipes_toolbox % remove the coll_v1, varcoll_v1, and po_v2 toolboxes from the search path
