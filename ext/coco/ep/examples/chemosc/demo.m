%% Bifurcation analysis of chemical oscillator
%
% We consider a model of oxidation of carbon monoxide on platinum analyzed
% in Tutorial IV: Two-parameter bifurcation analysis of equilibria and
% limit cycles with MATCONT, by Yu.A. Kuznetsov, September 20, 2011 (see
% http://www.staff.science.uu.nl/~kouzn101/NBA/LAB4.pdf). The analysis
% includes continuation of several one-dimensional branches of equilibria,
% as well as curves of saddle-node and Hopf bifurcations.

% Figure 1 shows several one-dimensional families of equilibria under
% variations in p2, for different values of p7, as well as projections onto
% the (p2,x1) plane of curves of saddle-node and Hopf bifurcations under
% simultaneous variations in p2 and p7. Figure 2 shows the corresponding
% bifurcation curves projected onto the (p2,p7) plane.

%% Initial encoding

% The continuation problem encoded below includes seven monitor functions
% that evaluate to each of the seven problem parameters, respectively, and
% seven corresponding inactive continuation parameters 'p1', 'p2' ,'p3',
% 'p4', 'p5', 'p6', and 'p7'. Its dimensional deficit equals 0. The call to
% the coco entry-point function indicates a desired manifold dimension of
% 1. To this end, the continuation parameter 'p2' is released and allowed
% to vary during continuation.

% Compute family of equilibria with p7=0.4

x0     = [0.001137; 0.891483; 0.062345];
pnames = {'p1' 'p2' 'p3' 'p4' 'p5' 'p6' 'p7'};
p0     = [2.5; 2.204678; 10; 0.0675; 1; 0.1; 0.4];

prob = coco_prob();
prob = coco_set(prob, 'ep', 'NSA', true); % Detect and locate neutral saddles

fprintf('\n Run=''%s'': Continue equilibria with p7=0.4.\n', 'p7=0.4');

bd1  = coco(prob, 'p7=0.4', 'ode', 'isol', 'ep', ...
  @bykov, @bykov_dx, @bykov_dp, x0, pnames, p0, ... % ep toolbox arguments
  1, 'p2', [0.4 3]);                                % cont toolbox arguments

% Repeat continuation with p7=0.15

p0(7) = 0.15;

fprintf('\n Run=''%s'': Continue equilibria with p7=0.15.\n', 'p7=0.15');

bd2 = coco(prob, 'p7=0.15', 'ode', 'isol', 'ep', ...
  @bykov, @bykov_dx, @bykov_dp, x0, pnames, p0, ... % ep toolbox arguments
  1, 'p2', [0.4 3]);                                % cont toolbox arguments

% Repeat continuation with p7=2.0

p0(7) = 2.0;

fprintf('\n Run=''%s'': Continue equilibria with p7=2.0.\n', 'p7=2.0');

bd3 = coco(prob, 'p7=2.0', 'ode', 'isol', 'ep', ...
  @bykov, @bykov_dx, @bykov_dp, x0, pnames, p0, ... % ep toolbox arguments
  1, 'p2', [0.4 3]);                                % cont toolbox arguments

%% Start continuation along curve of saddle-node bifurcations

% The continuation problem encoded below includes seven monitor functions
% that evaluate to each of the seven problem parameters, respectively, and
% seven corresponding inactive continuation parameters 'p1', 'p2' ,'p3',
% 'p4', 'p5', 'p6', and 'p7'. Its dimensional deficit equals -1. The call
% to the coco entry-point function indicates a desired manifold dimension
% of 1. To this end, the continuation parameters 'p2' and 'p7' are both
% released and allowed to vary during continuation.

labs = coco_bd_labs(bd1, 'SN');

fprintf(...
  '\n Run=''%s'': Continue SN points from point %d in run ''%s''.\n', ...
 'SN-curve', labs(2), 'p7=0.4');

bd4 = coco('SN-curve', 'ode', 'SN', 'SN', ...
  'p7=0.4', labs(2), ...    % ep toolbox arguments
  1, {'p2' 'p7'}, [0.4 3]); % cont toolbox arguments

%% Start continuation along curve of Hopf bifurcations

% The continuation problem encoded below includes seven monitor functions
% that evaluate to each of the seven problem parameters, respectively, and
% seven corresponding inactive continuation parameters 'p1', 'p2' ,'p3',
% 'p4', 'p5', 'p6', and 'p7'. Its dimensional deficit equals -1. The call
% to the coco entry-point function indicates a desired manifold dimension
% of 1. To this end, the continuation parameters 'p2' and 'p7' are both
% released and allowed to vary during continuation. In addition, a
% non-embedded monitor function and corresponding continuation parameter
% 'L1' are included and allow for the detection and location of generalized
% Hopf bifurcations corresponding to the transition between sub- and
% supercriticality.

labs = coco_bd_labs(bd1, 'HB');

prob = coco_prob();
prob = coco_set(prob, 'ep', 'BT', true); % detect and locate Bogdanov-Takens point
prob = coco_set(prob, 'cont', 'PtMX', 50);
prob = ode_HB2HB(prob, '', 'p7=0.4', labs(2));
[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
prob = coco_add_func(prob, 'lyap', @lyapunov, data.ep_eqn, ...
  'regular', 'L1', 'uidx', uidx);
prob = coco_add_event(prob, 'GH', 'L1', 0);

fprintf(...
  '\n Run=''%s'': Continue HB points from point %d in run ''%s''.\n', ...
 'HB-curve', labs(2), 'p7=0.4');

bd5 = coco(prob, 'HB-curve', [], ...  % fully constructed continuation problem
  1, {'p2' 'p7'}, [0.4 3]);           % cont toolbox arguments

%% Graphical representation of stored solutions

% One-parameter bifurcation diagram

figure(1); clf; hold on
thm = struct('special', {{'HB' 'SN' 'NSA'}});
coco_plot_bd(thm, 'p7=0.4', 'p2', 'x')
coco_plot_bd(thm, 'p7=0.15', 'p2', 'x')
coco_plot_bd(thm, 'p7=2.0', 'p2', 'x')
coco_plot_bd('SN-curve', 'p2', 'x')
thm.special = {'GH' 'BTP'};
thm.GH = {'kp', 'MarkerFaceColor', 'm', 'MarkerSize', 10};
coco_plot_bd(thm, 'HB-curve', 'p2', 'x')
hold off; grid on; axis([0.5 2 0 0.16])

% "Two-parameter bifurcation diagram
figure(2); clf; hold on
coco_plot_bd('SN-curve')
coco_plot_bd(thm, 'HB-curve')
hold off; grid on; axis([0.7 1.5 0.1 1.1])
