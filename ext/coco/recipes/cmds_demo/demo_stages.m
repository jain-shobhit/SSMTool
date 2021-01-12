%% Example 3.11

% The continuation problem encoded below in prob corresponds to a family of
% two zero functions in three continuation variables, a single monitor
% function, and a corresponding inactive continuation parameter 'p'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 0.

prob = coco_prob();
prob = coco_add_func(prob, 'fun1', @circ, [], 'zero', ...
  'u0', [0.9; 1.1]);
prob = coco_add_func(prob, 'fun2', @plan, [], 'zero', ...
  'uidx', 2, 'u0', -1.1);
prob = coco_add_func(prob, 'fun3', @hype, [], 'inactive', 'p', ...
  'uidx', [1; 3]);
coco(prob, 'run4', [], 0);

[data chart] = coco_read_solution('fun2', 'run4', 1);

%% Example 3.12

% The continuation problem encoded below in prob corresponds to a family of
% two zero functions in three continuation variables, a family of two
% monitor functions, an inactive continuation parameter 'p', and an active
% continuation parameter 'u3' tracking the value of the third continuation
% variable. Its dimensional deficit equals 0. The call to the coco
% entry-point function indicates a desired manifold dimension of 1. To this
% end, the continuation parameter 'p' is released and allowed to vary for
% the duration of the call to the coco entry-point function. Variations in
% 'u3' are output to screen (they are saved to the bd.mat file by
% default).

prob = coco_add_pars(prob, 'pars', 3, {'u3'}, 'active');
coco(prob, 'run5', [], 1, {'p' 'u3'}, [-0.4 0.4]);
