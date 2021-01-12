%% Example 15.1

% The continuation problem encoded below consists of a single zero function
% in terms of two continuation variables, a single monitor function, and
% the corresponding inactive continuation parameter 'p'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'p' is released and allowed to vary during continuation. Special points
% associated with the values 0.5, 1.0, 1.5, and 2.0 of 'p' are detected and
% identified with the event type 'UZ'.

prob = coco_add_func(coco_prob(), 'circ', @circ, [], 'zero', ...
  'u0', [1; 1.1]);
prob = coco_add_func(prob, 'fun2', @euclid, [], 'inactive', 'p', ...
  'uidx', [1; 2]); % [In Recipes for Continuation, function handle was @dist]
prob = coco_add_event(prob, 'UZ', 'p', 0.5:0.5:2);
coco(prob, 'run', [], 1, 'p', [0.5, 3]);

%% Example 15.1 continued [Not in Recipes for Continuation]

% The continuation problem structure encoded below has an additional two
% monitor functions that evaluate to the continuation variables, and the
% corresponding active continuation parameters 'u1' and 'u2'. Its
% dimensional deficit still equals 0. A one-dimensional family of solutions
% is obtained by again releasing 'p' and allowing it to vary during
% continuation. The additional continuation parameters are stored with the
% bifurcation data, thereby allowing straightforward extraction and
% plotting.

prob = coco_add_pars(prob, 'vars', [1, 2], {'u1', 'u2'}, 'active');
bd = coco(prob, 'run', [], 1, {'p', 'u1', 'u2'}, [0.1, 5]);

figure(1)
clf
hold on
th=0:0.01:2*pi;
plot(cos(th),1+sin(th),'g')

x   = coco_bd_col(bd, 'u1');  % Extract column data
y   = coco_bd_col(bd, 'u2');  % Extract column data
idx = coco_bd_idxs(bd, 'UZ'); % Extract row indices
plot(x(idx), y(idx), 'go', 'LineWidth', 2, 'MarkerSize', 10);

for r = 0.5:0.5:2
    thmin=atan2(r/2,sqrt(4-r^2)/2);
    thmax=atan2(r/2,-sqrt(4-r^2)/2);
    th=thmin:0.01:thmax;
    plot(r*cos(th),r*sin(th),'k-', 'LineWidth', 2)
end
hold off
axis equal
