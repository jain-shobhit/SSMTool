coco_use_recipes_toolbox atlas1d_v8 % Add atlas1d_v8 atlas algorithm to search path

%% Example 20.4 

% The continuation problem structure encoded below consists of N zero
% functions in terms of N+1 continuation variables (N mesh values and one
% problem parameter), a single monitor function that evaluates the problem
% parameter, and the corresponding inactive continuation parameter 'p'. Its
% dimensional deficit equals 0. A one-parameter family of soutions is
% obtained by releasing 'p' and allowing it to vary during continuation.
% Solutions correspond to the values of the function tanh(p*t)/tanh(p) on a
% uniform mesh without adaptation to the local properties of the function.

N          = 7;
p0         = 1;
data       = struct();
data.x_idx = 1:N;
data.p_idx = N+1;
data.t     = linspace(-1, 1, N)'; % Uniform mesh
u0         = [tanh(p0*data.t)/tanh(p0); p0];

prob = coco_add_func(coco_prob(), 'tanh', @tanh_F, data, ...
  'zero', 'u0', u0);
prob = coco_add_slot(prob, 'tanh', @coco_save_data, data, 'save_full');
prob = coco_add_pars(prob, 'pars', N+1, 'p');
prob = coco_add_event(prob, 'UZ', 'p', [1 4 7 10]);
prob = coco_set(prob, 'cont', 'h0', 1', 'hmax', 1);
coco(prob, 'run1', [], 1, 'p', [1 11]);

%% Example 20.5

% The continuation problem structure encoded below consists of N zero
% functions in terms of N+1 continuation variables (N mesh values and one
% problem parameter), a single monitor function that evaluates the problem
% parameter, and the corresponding inactive continuation parameter 'p'. Its
% dimensional deficit equals 0. A one-parameter family of soutions is
% obtained by releasing 'p' and allowing it to vary during continuation.
% Solutions correspond to the values of the function tanh(p*t)/tanh(p) on a
% moving mesh adapted to the local properties of the function. For the
% given adaptation window, the discretization order remains constant.

data.N    = N;
data.pdim = 1;
data.xtr  = zeros(N+data.pdim,1);
data.xtr([1 N:N+data.pdim]) = [1 N:N+data.pdim]; % Index array of invariant elements
data.th   = linspace(-1, 1, N)'; % Equidistributed reference mesh
data.s    = .3; % Scaling
data.HINC = 2;  % Upper bound on adaptation window
data.HDEC = 0;  % Lower bound on adaptation window
data = coco_func_data(data); % Convert to func_data for shared access

data.t    = [-1; 1/p0*atanh((-1+2/(N-1)*(1:N-2)')*tanh(p0)); 1]; % Initial equidistributed mesh
u0        = [tanh(p0*data.t)/tanh(p0); p0];
prob = coco_add_func(coco_prob(), 'tanh', @tanh_F, data, 'zero', ...
  'u0', u0, 'remesh', @remesh); % Include 'remesh' option
prob = coco_add_slot(prob, 'tanh', @coco_save_data, data, 'save_full');
prob = coco_add_pars(prob, 'pars', N+1, 'p');
prob = coco_add_event(prob, 'UZ', 'p', [1 4 7 10]);
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'h0', 1, 'hmax', 1); % Remesh after each continuation step
coco(prob, 'run2', [], 1, 'p', [1 11]);

%% Example 20.5

% The continuation problem structure encoded below consists of N zero
% functions in terms of N+1 continuation variables (N mesh values and one
% problem parameter), a single monitor function that evaluates the problem
% parameter, and the corresponding inactive continuation parameter 'p'. Its
% dimensional deficit equals 0. A one-parameter family of soutions is
% obtained by releasing 'p' and allowing it to vary during continuation.
% Solutions correspond to the values of the function tanh(p*t)/tanh(p) on a
% moving mesh adapted to the local properties of the function. For the
% given adaptation window, the discretization order varies during
% continuation.

data.t    = [-1; 1/p0*atanh((-1+2/(N-1)*(1:N-2)')*tanh(p0)); 1]; % Initial equidistributed mesh
u0        = [tanh(p0*data.t)/tanh(p0); p0];
data.N    = N;
data.pdim = 1;
data.xtr  = zeros(N+data.pdim,1);
data.xtr([1 N:N+data.pdim]) = [1 N:N+data.pdim]; % Index array of invariant elements
data.th   = linspace(-1, 1, N)'; % Equidistributed reference mesh
data.s    = 0.3; % Scaling
data.HINC = 0.3; % Upper bound on adaptation window
data.HDEC = 0.2; % Lower bound on adaptation window
data = coco_func_data(data); % Convert to func_data for shared access

prob = coco_add_func(coco_prob(), 'tanh', @tanh_F, data, 'zero', ...
  'u0', u0, 'remesh', @remesh); % Include 'remesh' option
prob = coco_add_slot(prob, 'tanh', @coco_save_data, data, 'save_full');
prob = coco_add_pars(prob, 'pars', N+1, 'p');
prob = coco_add_event(prob, 'UZ', 'p', [1 4 7 10]);
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'h0', 1, 'hmax', 1); % Remesh after each continuation step
coco(prob, 'run3', [], 1, 'p', [1 11]);

coco_use_recipes_toolbox % Remove atlas1d_v8 atlas algorithm from search path
