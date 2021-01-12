%% Example 3.13, first version of period

% The continuation problem encoded below in prob corresponds to a family of
% eight zero functions in ten continuation variables, a family of two
% monitor functions, and the corresponding inactive continuation parameters
% 'a' and 'b'. Its dimensional deficit equals 0. The call to the coco
% entry-point function indicates a desired manifold dimension of 1. To this
% end, the continuation parameter 'a' is released and allowed to vary
% during continuation.

u0 = [1; 0.3; 1.275; -0.031; -0.656; 0.382; 0.952; ...
  -0.197; -0.103; 0.286];
prob = period_A(u0, 4);
prob = coco_add_pars(prob, 'pars', [1 2], {'a' 'b'});
coco(prob, 'run6A', [], 1, 'a', [0.8 1.2]);

%% Example 3.13, second version of period

% The continuation problem encoded below in prob corresponds to a family of
% ten zero functions in twelve continuation variables, a family of two
% monitor functions, and the corresponding inactive continuation parameters
% 'a' and 'b'. Its dimensional deficit equals 0. The call to the coco
% entry-point function indicates a desired manifold dimension of 1. To this
% end, the continuation parameter 'a' is released and allowed to vary
% during continuation.

u0 = [1; 0.3; 1.275; -0.031; -0.656; 0.382; 0.952; ...
  -0.197; -0.103; 0.286; 1.275; -0.031];
prob = period_B(u0, 4);
prob = coco_add_pars(prob, 'pars', [1 2], {'a' 'b'});
coco(prob, 'run6B', [], 1, 'a', [0.8 1.2]);
