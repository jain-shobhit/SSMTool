%% A web of bifurcation curves for a three-dimensional vector field
%
% We explore bifurcations of equilibria and periodic orbits in a vector
% field described in Freire, E., Rodriguez-Luis, A., Gamero, E. & Ponce, E.
% (1993), ''A case study for homoclinic chaos in an autonomous electronic
% circuit: A trip from Takens-Bogdanov to Hopf-Shilnikov,'' Physica D 62,
% 230?253, and also used in the AUTO manual to demonstrate functionality
% and syntax.

% Figure shows trival branch of equilibria at the origin under variations
% in one problem parameter, the detection and continuation under variations
% in two problem parameters of a family of Hopf bifurcations of the trivial
% equilibrium, a primary branch of periodic orbits emanating from the Hopf
% bifurcation under variations in one problem parameter, a secondary branch
% of periodic orbits emanating from a branch point along the primary
% branch, continuation of families of saddle-node, period-doubling and
% torus bifurcations under variations in two problem parameters from points
% along the secondary branf of periodic orbits, and continuation of a
% branch of period-doubled periodic orbits under variation in one problem
% parameter emanating from a point on the period-doubling bifurcation
% branch.

%% Continue branch of equilibria at origin

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'nu' is released and allowed to vary during
% continuation.

%       nu    be     ga     r     a3    b3    
p0 = [-0.65 ; 0.5 ; -0.6 ; 0.6 ; 0.3 ; 0.9];

prob = coco_prob();
prob = ode_isol2ep(prob, '', @tor, [0;0;0], ...
  {'nu', 'be', 'ga', 'r', 'a3', 'b3'}, p0);

fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep');

bd_ep = coco(prob, 'ep', [], 1, 'nu', [-0.65, -0.55]);

%% Continue branch of Hopf bifurcations of equilibrium at origin

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'nu' and 'be' are released and allowed to vary
% during continuation.

lab = coco_bd_labs(bd_ep, 'HB');
prob = coco_prob();
prob = ode_HB2HB(prob, '', 'ep', lab);

fprintf(...
  '\n Run=''%s'': Continue Hopf bifurcations from point %d in run ''%s''.\n', ...
  'hb', lab, 'ep');

bd_hb = coco(prob, 'hb', [], 1, {'nu', 'be'}, [-0.65, -0.55]);

%% Continue first branch of periodic orbits

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'nu' is released and allowed to vary during
% continuation. The screen output also includes the orbital period, the
% ratio of the interpolation error estimate with the tolerance, and the
% stability indicator that counts the number of Floquet multipliers outside
% the unit circle.

prob = coco_prob();
prob = ode_HB2po(prob, '', 'ep', lab);
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'PtMX', [100 0]);

fprintf(...
  '\n Run=''%s'': Continue primary branch of periodic orbits from point %d in run ''%s''.\n', ...
  'po1', lab, 'ep');

bd1  = coco(prob, 'po1', [], 1, ...
  {'nu' 'po.period' 'po.orb.coll.err_TF' 'po.test.USTAB'}, [-0.65, -0.55]);

%% Switch branch at BP

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'nu' is released and allowed to vary during
% continuation. The screen output also includes the orbital period, the
% ratio of the interpolation error estimate with the tolerance, and the
% stability indicator that counts the number of Floquet multipliers outside
% the unit circle.

BPlabs = coco_bd_labs(bd1, 'BP');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'PtMX', [100 0]);

fprintf(...
  '\n Run=''%s'': Continue secondary branch of periodic orbits from point %d in run ''%s''.\n', ...
  'po2', BPlabs(end), 'po1');

bd2 = coco(prob, 'po2', 'ode', 'BP', 'po', 'po1', BPlabs(1), 1, ...
  {'nu' 'po.period' 'po.orb.coll.err_TF' 'po.test.USTAB'}, [-0.65, -0.55]);

%% Continuation of saddle-node bifurcations

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'nu' and 'be' are released ('po.period' is
% already active) and allowed to vary during continuation.

SNlabs = coco_bd_labs(bd1, 'SN');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5);

fprintf(...
  '\n Run=''%s'': Continue saddle-node bifurcations from point %d in run ''%s''.\n', ...
  'lp1', SNlabs(2), 'po1');

bd3 = coco(prob, 'lp1', 'ode', 'SN', 'SN', ...
	'po1', SNlabs(2), {'nu' 'po.period' 'be'}, [-0.65, -0.55]);

%% Continuation of period-doubling bifurcations

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'nu' and 'be' are released ('po.period' is
% already active) and allowed to vary during continuation.

PDlabs = coco_bd_labs(bd2, 'PD');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5);

fprintf(...
  '\n Run=''%s'': Continue period-doubling bifurcations from point %d in run ''%s''.\n', ...
  'pd1', PDlabs(1), 'po2');

bd4 = coco(prob, 'pd1', 'ode', 'PD', 'PD', ...
  'po2', PDlabs(1), {'nu' 'po.period' 'be'}, [-0.65, -0.55]);

%% Continuation of torus bifurcations

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'nu' and 'be' are released ('po.period' is
% already active) and allowed to vary during continuation.

TRlabs = coco_bd_labs(bd2, 'TR');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5);

fprintf(...
  '\n Run=''%s'': Continue torus bifurcations from point %d in run ''%s''.\n', ...
  'tr1', TRlabs(1), 'po2');

bd5 = coco(prob, 'tr1', 'ode', 'TR', 'TR', ...
  'po2', TRlabs(1), {'nu' 'po.period' 'be'}, [-0.65, -0.55]);

%% Continuation along period-doubled branch

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'nu' is released and allowed to vary during
% continuation. The screen output also includes the orbital period and the
% ratio of the interpolation error estimate with the tolerance.

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5);

fprintf(...
  '\n Run=''%s'': Continue period-doubled periodic orbits from point %d in run ''%s''.\n', ...
  'po_db', 4, 'pd1');

bd_db  = coco(prob, 'po_db', 'ode', 'PD', 'po', 'pd1', 4, 1, ...
  {'nu' 'po.period' 'po.orb.coll.err_TF'}, [-0.65, -0.55]);

%% Graphical representation of results

figure(1); clf; hold on; grid on
thm = struct('special', {{'EP', 'HB'}});
coco_plot_bd(thm, 'ep', 'nu', 'be', '||x||_2')
thm = struct('special', {{'EP', 'BP', 'FP'}});
coco_plot_bd(thm, 'hb', 'nu', 'be', '||x||_2')
thm = struct('special', {{'EP', 'SN', 'BP', 'FP'}});
coco_plot_bd(thm, 'po1', 'nu', 'be', '||x||_{2,MPD}')
thm = struct('special', {{'EP', 'PD', 'TR'}});
coco_plot_bd(thm, 'po2', 'nu', 'be', '||x||_{2,MPD}')
thm = struct('special', {{'EP', 'FP'}});
coco_plot_bd(thm, 'lp1', 'nu', 'be', '||x||_{2,MPD}')
thm = struct('special', {{'EP', 'FP'}});
coco_plot_bd(thm, 'pd1', 'nu', 'be', '||x||_{2,MPD}')
thm = struct('special', {{'EP', 'FP'}});
coco_plot_bd(thm, 'tr1', 'nu', 'be', '||x||_{2,MPD}')
thm = struct('special', {{'EP', 'BP', 'FP', 'PD'}}, 'xlab', '\nu', 'ylab', '\beta');
coco_plot_bd(thm, 'po_db', 'nu', 'be', '||x||_{2,MPD}')
hold off; view(-192,44)
