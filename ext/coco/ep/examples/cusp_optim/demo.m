%% Stationary points in the cusp normal form
%
% This demo illustrates the successive continuation approach to searching
% for stationary points that satisfy the necessary conditions for a local
% extremum in the variable kappa along the two-dimensional manifold of
% equilibria of the cusp normal form x' = kappa-x*(lambda-x^2).
%
% The augmented continuation problem constructed below has dimensional
% deficit 0-2=-2. The continuation parameters 'ka', 'la', 'd.ka', and
% 'd.la' are all initially inactive. Since each stage of continuation is
% along a one-dimensional manifold, three of these parameters must be
% released in each stage. In particular, 'ka' must always be active.
% Consequently, 'd.ka' should be active during the first two stages of
% continuation, the second of which terminates when this parameter equals
% 1. It should be inactive (and equal to 1) for all remaining stages of
% continuation. Since 'la' is inactive during the first two stages of
% continuation and released only in the third stage of continuation, 'd.la'
% should be active during all three stages.
%
% In the first stage of continuation, local extrema in 'ka' are located
% along a one-dimensional solution manifold with trivial Lagrange
% multipliers. Each of these is a branch point, from which emanates a
% secondary one-dimensional submanifold along which the Lagrange
% multipliers take on nontrivial values. As explained above, we terminate
% continuation along such a manifold when 'd.ka' equals 1. Stationary
% points within the computational domain correspond to points with
% vanishing 'd.la'.
%
% In this example, the third stage of continuation is along a family of
% saddle-node bifurcations with the cusp coincident with the stationary
% point.

% The figure shows the one-dimensional solution manifolds obtained in the
% first and third stages of continuation.

%% Initial encoding

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', false);

fcns = {@cusp, @cusp_dx, @cusp_dp, @cusp_dxdx, @cusp_dxdp, @cusp_dpdp};

%% First run to find local extrema
% zero problem
prob1 = ode_isol2ep(prob, '', fcns{:}, 0, {'ka' 'la'}, [0; 0.5]);

% adjoint
prob1 = adjt_isol2ep(prob1, '');

% continuation
cont_args = {1, {'ka' 'd.ka' 'd.la'}, [-0.5 0.5]};
bd1 = coco(prob1, 'cusp1', [], cont_args{:});

%% Switch at fold to branch with nontrivial multipliers
BPlab = coco_bd_labs(bd1, 'BP');

% zero problem
prob2 = ode_BP2ep(prob, '', 'cusp1', BPlab(1));

% adjoint
prob2 = adjt_BP2ep(prob2, '', 'cusp1', BPlab(1));

% continuation
cont_args = {1, {'d.ka' 'ka' 'd.la'}, {[0 1] [-0.5 0.5]}};
bd2 = coco(prob2, 'cusp2', [], cont_args{:});

%% Constrain nontrivial multiplier and release additional continuation parameters
lab = coco_bd_labs(bd2, 'EP');

% zero problem
prob3 = ode_ep2ep(prob, '', 'cusp2', lab(2));

% adjoint
prob3 = adjt_ep2ep(prob3, '', 'cusp2', lab(2));

% events
prob3 = coco_add_event(prob3, 'OPT', 'd.la', 0);

% continuation
cont_args = {1, {'d.la' 'ka' 'la'}, {[], [-0.5 0.5], [-2 2]}};
coco(prob3, 'cusp3', [], cont_args{:});

%% Graphical representation

figure(1); clf; hold on
coco_plot_bd('cusp1', 'ka', 'la', 'x')
thm         = struct();
thm.ustab   = '';
thm.lspec   = {'g-', 'LineWidth', 1};
thm.special = {'OPT'};
thm.OPT     = {'kp', 'MarkerFaceColor', 'r', 'MarkerSize', 8};
coco_plot_bd(thm, 'cusp3', 'ka', 'la', 'x')
hold off; grid on; view(3)
