%% Example 3.7

% The continuation problem encoded below in prob corresponds to a single
% zero function in two continuation variables. Its dimensional deficit
% equals 1. The call to the coco entry-point function indicates a desired
% manifold dimension of 1.

prob = coco_prob();
prob = coco_add_func(prob, 'fun1', @circ, [], 'zero', ...
  'u0', [0.9; 1.1]);
coco(prob, 'run1', [], 1);

%% Example 3.9

% The continuation problem encoded below in prob corresponds to a single
% zero function in two continuation variables, a single monitor function,
% and a corresponding inactive continuation parameter 'p'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 0.

prob = coco_add_func(prob, 'fun2', @dist, [], 'inactive', 'p', ...
  'uidx', [1; 2]);
coco(prob, 'run2', [], 0);

%% Example 3.10

% The continuation problem encoded in prob still corresponds to a single
% zero function in two continuation variables, a single monitor function,
% and a corresponding inactive continuation parameter. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'p' is released and allowed to vary during continuation.

coco(prob, 'run3', [], 1, 'p', [0.1 5]);

%% Reading from disk

% The coco_read_solution utility extracts the chart structure and data
% array associated with a given function identifier from a solution file.
% The coco_bd_read utility extracts bifurcation data from the bd.mat file
% in the corresponding data directory.

[data chart] = coco_read_solution('fun1', 'run3', 7); %Extract solution from disk
chart.x

bd = coco_bd_read('run3'); % Extract bifurcation data from disk
bd(1:8,:)
