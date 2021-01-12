%% Example 3.14

% The continuation problem encoded below in prob corresponds to a family of
% 2N zero functions in 2N+6 continuation variables, a family of 6 monitor
% functions, and the corresponding inactive continuation parameters 'f0',
% 'g0', 'f1', 'g1', 'p1', and 'p2'. Its dimensional deficit equals 0. The
% call to the coco entry-point function indicates a desired manifold
% dimension of 1. To this end, the continuation parameter 'g0' is released
% and allowed to vary during continuation. The assignment of the variable
% data as function data structure for the coco_save_data slot function
% ensures that the content of data is stored to the data array during
% continuation.

N = 40;
data.dep_idx = (1:2*N+4)';
data.par_idx = [2*N+5; 2*N+6];
data.f_idx   = (1:2:2*N+3)';
data.g_idx   = (2:2:2*N+4)';
oneN   = ones(2*N,1);
zeroN  = zeros(2*N,1);
data.A = spdiags([oneN zeroN -2*oneN zeroN oneN], ...
  [0 1 2 3 4], 2*N, 2*N+4);
data.B = 1/(12*(N+1)^2)*spdiags([oneN zeroN 10*oneN zeroN oneN], ...
  [0 1 2 3 4], 2*N, 2*N+4);
x0   = ones(1,N+2);
y0   = ones(1,N+2);
p0   = [1; 1];
dep0 = [x0; y0];
u0   = [dep0(:); p0];
prob = coco_prob();
prob = coco_add_func(prob, 'finitediff', @finitediff, data, ...
  'zero', 'u0', u0);

prob = coco_add_pars(prob, 'pars', [1:2 2*N+3:2*N+6]', ...
  {'f0' 'g0' 'f1' 'g1' 'p1' 'p2'});

% The following line is not in Recipes for Continuation, but included here to support
% graphing in figure_3_3.m
prob = coco_add_slot(prob, 'finitediff', @coco_save_data, data, 'save_full');

coco(prob, 'brusselator', [], 1, 'g0', [0 10]);

%% Example 4.9

% The continuation problem encoded in prob includes a zero or monitor
% function object with function identifier 'finitediff'. The call below to
% the coco_get_func_data utility returns the corresponding function data
% structure.

coco_get_func_data(prob, 'finitediff', 'data')
