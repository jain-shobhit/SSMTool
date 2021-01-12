%% The cusp normal form from Section 16.2.3 of Recipes for Continuation
%
% The cusp normal form is given by the one-dimensional, autonomous vector
% field f = kappa-x*(lambda-x^2) in terms of the state variable x and the
% problem parameters kappa and lambda. We rely on an anonymous,
% non-vectorized encoding of the vector field and perform several stages of
% continuation of equilibria and saddle-node bifurcations.

% Figure 1 shows a one-dimensional family of equilibria under variations
% in kappa, for fixed lambda. Saddle-node bifurcations are identified by
% green circles. Figure 2 shows a one-dimensional family of saddle-node
% equilibria under simultaneous variation in kappa and lambda. The figure
% shows the result of starting continuation from a saddle-node bifurcation
% detected during equilibrium continuation, as well as from a solution to
% the saddle-node equilibrium continuation problem.

%% Initial encoding

% The continuation problem encoded below includes two monitor functions
% that evaluate to the problem parameters kappa and lambda, respectively,
% and two corresponding inactive continuation parameters 'ka' and 'la'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'ka' is released and allowed to vary during
% continuation.

% Define cusp normal form using inline, anonymous functions.

cusp    = @(x,p) p(1)-x*(p(2)-x^2);
cusp_dx = @(x,p) 3*x^2-p(2);
cusp_dp = @(x,p) [1 -x];

% Initialize continuation problem and settings associated with toolbox
% constructor.

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', false);

% Compute family of equilibria using nested call to ode_isol2ep. Note:
% derivatives can be omitted, see help for ode_isol2ep.

fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  'cusp1');

bd1 = coco(prob, 'cusp1', 'ode', 'isol', 'ep', ...
  cusp, cusp_dx, cusp_dp, 0, {'ka' 'la'}, [0; 0.5], ... % ep toolbox arguments
  1, {'ka' 'la'}, [-0.5 0.5]);                          % cont toolbox arguments

%% Start continuation along solutions to saddle-node equilibrium problem

% The continuation problem encoded below includes two monitor functions
% that evaluate to the problem parameters kappa and lambda, respectively,
% and two corresponding inactive continuation parameters 'ka' and 'la'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'ka' and 'la' are both released and allowed to
% vary during continuation.

% Extract solution label associated with saddle-node bifurcations.

labs = coco_bd_labs(bd1, 'SN');

% Compute family of saddle-node equilibria using nested call to ode_SN2SN.
% Note: an equivalent construction results from using the ode_ep2SN
% constructor. See help for ode_SN2SN.

fprintf(...
  '\n Run=''%s'': Continue SN points from point %d in run ''%s''.\n', ...
  'cusp2', labs(1), 'cusp1');

bd2 = coco(prob, 'cusp2', 'ode', 'SN', 'SN', ...
  'cusp1', labs(1),      ...   % ep toolbox arguments
  1, {'ka' 'la'}, [-0.5 0.5]); % cont toolbox arguments

%% Test restart of saddle-node continuation

% The continuation problem encoded below includes two monitor functions
% that evaluate to the problem parameters kappa and lambda, respectively,
% and two corresponding inactive continuation parameters 'ka' and 'la'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'ka' and 'la' are both released and allowed to
% vary during continuation.

fprintf(...
  '\n Run=''%s'': Continue SN points from point %d in run ''%s''.\n', ...
  'cusp3', 2, 'cusp2');

bd3 = coco(prob, 'cusp3', 'ode', 'SN', 'SN', ...
  'cusp2', 2,            ...  % ep toolbox arguments
  1, {'ka' 'la'}, [-1 -0.5]); % cont toolbox arguments

%% Graphical representation of stored solutions

% "One-parameter continuation" of equilibria
figure(1); clf
thm = struct('special', {{'SN'}});
coco_plot_bd(thm, 'cusp1', 'ka', 'x')
grid on

% "Two-parameter continuation" of saddle-node bifurcations
figure(2); clf; hold on
thm = struct('special', {{'EP'}});
coco_plot_bd(thm, 'cusp2')
coco_plot_bd(thm, 'cusp3')
hold off; grid on
