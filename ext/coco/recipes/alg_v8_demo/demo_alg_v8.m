coco_use_recipes_toolbox alg_v8 % add the alg_v8 toolbox to the search path

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
bd1 = coco(prob, 'regular', [], 1, {'ka' 'la' 'alg.test.FO'}, [-0.5 0.5]);

figure(1)
clf
hold on

subplot(1,2,1)
x  = coco_bd_col(bd1, '||U||'); % Extract column data
ka = coco_bd_col(bd1, 'ka'); % Extract column data
plot(ka, x, '.-')
hold on
idx = coco_bd_idxs(bd1, 'EP'); % Extract row indices
plot(ka(idx),x(idx),'go');
idx = coco_bd_idxs(bd1, 'FO'); % Extract row indices
plot(ka(idx),x(idx),'ko');
hold off
grid on
drawnow

%% Continuation of fold points with Moore-Spence system

% The continuation problem encoded below consists of three zero functions
% (the 'alg' zero function and the Moore-Spence system) in terms of four
% continuation variables (the problem variables, problem parameters, and
% nullvector), two embedded monitor functions, and the corresponding inactive
% continuation parameters 'ka' and 'la'. Its dimensional deficit equals -1.
% The call to the coco entry-point function indicates a desired manifold
% dimensionality of 1. To this end, the 'ka' and 'la' continuation
% parameters are released and allowed to vary during continuation.

labs = coco_bd_labs(bd1, 'FO'); % Extract solution label
prob = coco_prob();
prob = alg_FO2FO(prob, '', 'regular', labs(1));
bd2  = coco(prob, 'cusp', [], 1, {'ka' 'la'}, {[-2 2] [-2 2]});

% plot bifurcation diagrams
subplot(1,2,2)
ka = coco_bd_col(bd2, 'ka'); % Extract column data
la = coco_bd_col(bd2, 'la'); % Extract column data
plot(ka, la, '.-')
hold on
idx = coco_bd_idxs(bd2, 'EP'); % Extract row indices
plot(ka(idx),la(idx),'go');
hold off
grid on
drawnow

coco_use_recipes_toolbox % remove the alg_v8 toolbox from the search path
