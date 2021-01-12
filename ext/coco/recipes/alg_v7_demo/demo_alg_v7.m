coco_use_recipes_toolbox alg_v7 % add the alg_v7 toolbox to the search path

%% Run with nonembedded fold monitor function

% The continuation problem encoded below consists of a single zero function
% in terms of three continuation variables, two embedded monitor functions
% with corresponding inactive continuation parameters 'ka' and 'la', and a
% single nonembedded monitor function with corresponding continuation
% parameter 'alg.test.FO'. Its dimensional deficit equals 0. The call to
% the coco entry-point function indicates a desired manifold dimension of
% 1. To this end, the continuation parameter 'ka' is released and allowed
% to vary during continuation. Zero crossings of 'alg.test.FO' along the
% solution manifold corresponding to event type 'FO' are detected and
% located during continuation.

alg_args = {@cusp, @cusp_DFDX, @cusp_DFDP, 0, {'ka' 'la'}, ...
  [0; 0.5]};
prob = alg_isol2eqn(coco_prob(), '', alg_args{:});
coco(prob, 'regular', [], 1, {'ka' 'la' 'alg.test.FO'}, [-0.5 0.5]);

%% Run with embedded fold monitor function

% The continuation problem encoded below is identical to that above, but
% the monitor function associated with the active 'alg.test.FO'
% continuation parameter is embedded in the extended continuation problem.

alg_args = {@cusp, @cusp_DFDX, @cusp_DFDP, 0, {'ka' 'la'}, ...
  [0; 0.5]};
prob = coco_set(coco_prob(), 'alg', 'FO', 'active');
prob = alg_isol2eqn(prob, '', alg_args{:});
coco(prob, 'active', [], 1, {'ka' 'la' 'alg.test.FO'}, [-0.5 0.5]);

%% Restarting from fold point with embedded monitor function

% The continuation problem encoded below is identical to the second one
% above, but reconstructed from a fold point detected during the first run,
% and with 'la' active and 'alg.test.FO' inactive. 

bd1 = coco_bd_read('regular');
prob = coco_set(coco_prob(), 'alg', 'FO', 'active');
labs = coco_bd_labs(bd1, 'FO'); % Extract solution label
prob = alg_sol2eqn(prob, '', 'regular', labs(1));
prob = coco_xchg_pars(prob, 'la', 'alg.test.FO'); % Swith between II and JJ
coco(prob, 'cusp', [], 1, {'ka' 'la' 'alg.test.FO'}, [-0.5 0.5]);

%% Two-dimensional atlas algorithm with nonembedded fold monitor function
coco_use_recipes_toolbox atlas2d_v6 alg_v7 % add the atlas2d_v6 toolbox to the search path

% The continuation problem encoded below is identical to the first one
% above. The call to the coco entry-point function indicates a desired
% manifold dimension of 2. To this end, the continuation parameters 'ka'
% and 'la' are released and allowed to vary during continuation.

alg_args = {@cusp, @cusp_DFDX, @cusp_DFDP, 0, {'ka' 'la'}, ...
  [0; 0.5]};
prob = coco_prob();
prob = coco_set(prob, 'alg', 'norm', true); % Store the absolute value of x to the bifurcation data
prob = alg_isol2eqn(prob, '', alg_args{:});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create); % Identify the atlas algorithm
prob = coco_set(prob, 'cont', 'h', .075, 'almax', 35);
prob = coco_set(prob, 'cont', 'NPR', 100, 'PtMX', 2000);
coco(prob, 'cuspsurface', [], 2, {'ka' 'la' 'alg.test.FO'}, ...
  {[-0.5 0.5], [-1 1]});

coco_use_recipes_toolbox % remove the atlas2d_v6 and alg_v7 toolboxes from the search path
