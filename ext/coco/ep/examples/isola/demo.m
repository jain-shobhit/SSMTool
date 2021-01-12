%% Continuation of isolas of equilibrium points
%
% We consider a continuous stirred tank reactor model of chemical reactions
% analyzed in "On the Numerical Continuation of Isolas of Equilibria," by
% Avitabile, Desroches, and Rodriquez, International Journal of Bifurcation
% and Chaos, 22(11), art. no. 1250277, 2012. The analysis uncovers closed
% one-dimensional families of equilibria--known as isolas--and implements a
% continuation problem corresponding to the simultaneous continuation of a
% discretization of such isolas.

% Figure 1 shows several isolas of equilibria under simultaneous variations
% in the problem parameters, as well as a projection of a curves of Hopf
% bifurcations onto the (tau,x1) plane. Figure 2 shows the corresponding
% variations in weighted arclength.

%% Initial encoding

% The continuation problem encoded below includes two monitor functions
% that evaluate to each of the two problem parameters, respectively, and
% two corresponding inactive continuation parameters 'tau' and 'lambda'.
% Its dimensional deficit equals 0. The call to the coco entry-point
% function indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'tau' is released and allowed to vary during
% continuation.

% Compute initial isola family of equilibria

x0     = [0.75 4];
pnames = {'tau' 'lambda'};
p0     = [0.5; 0.11];

fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  'initial');

bd1  = coco('initial', 'ode', 'isol', 'ep', ... 
  @cstr, @cstr_dx, @cstr_dp, x0, pnames, p0, ... % ep toolbox arguments
  1, 'tau', [0 2]);                              % cont toolbox arguments

%% Start continuation of isola discretization, anchored at Hopf bifurcation

% The continuation problem encoded below includes three monitor functions
% that evaluate to each of the two problem parameters and the total
% weighted arclength, respectively, and three corresponding inactive
% continuation parameters 'tau', 'lambda', and 'L'. Its dimensional deficit
% equals -2. The call to the coco entry-point function indicates a desired
% manifold dimension of 1. To this end, the continuation parameters 'L',
% 'lambda', and 'tau' are all released and allowed to vary during
% continuation.

idxs = coco_bd_idxs(bd1, 'HB');
vars = coco_bd_col(bd1, 'x');
pars = coco_bd_col(bd1, {'tau' 'lambda'});

weights = [10; 1; 6.5; 1];

% Compute a weighted arclength along isola
s = 0;
for idx = idxs(2)+1:idxs(3)
dw = ( [vars(:,idx); pars(:, idx)] ...
	  - [vars(:,idx-1); pars(:, idx-1)] ).*weights;
s = [s; s(end) + norm(dw)]; %#ok<AGROW>
end

% Construct Hopf bifurcation continuation problem
labs = coco_bd_labs(bd1, 'HB');

prob = coco_prob();
prob = ode_HB2HB(prob, 'isola1', 'initial', '', labs(2));

% Construct zero problems for an interpolation of equilibria along the
% initial isola
N = 50;
vars = interp1(s, vars(:,idxs(2):idxs(3))', 0:s(end)/N:s(end))';
pars = interp1(s, pars(:,idxs(2):idxs(3))', 0:s(end)/N:s(end))';
prob = coco_set(prob, 'ep', 'SN', 'off', 'HB', 'off');
for idx = 2:N
  x0   = vars(:,idx);
  p0   = pars(:,idx);
  oid  = sprintf('isola%d', idx);
  prob = ode_isol2ep(prob, oid, @cstr, @cstr_dx, @cstr_dp, x0, p0);
end

% Glue parameter values and impose distance conditions along the isola
[data, uidx] = coco_get_func_data(prob, 'isola1.ep', 'data', 'uidx');
varidx = uidx(data.ep_eqn.x_idx);
paridx = uidx(data.ep_eqn.p_idx);
for idx = 2:N
  fid    = sprintf('isola%d.ep', idx);
  [data, uidx] = coco_get_func_data(prob, fid, 'data', 'uidx');
  varidx = [varidx uidx(data.ep_eqn.x_idx)]; %#ok<AGROW>
  paridx = [paridx uidx(data.ep_eqn.p_idx)]; %#ok<AGROW>
end

prob = coco_add_glue(prob, 'glue', paridx(2,1:end-1), paridx(2,2:end));
uidx = [varidx; paridx(1,:)];
data = struct('w', repmat(weights(1:3), [1 N]));
np   = numel(uidx)/3;
data.np  = np;
data.shp = [3 np];
data.shf = circshift(reshape(1:3*np, [3 np]), [0 -1]);
prob = coco_add_func(prob, 'dist', @wdist, data, ...
  'zero', 'uidx', uidx, 'u0', s(end)/N);
			   
% Introduce a continuation parameter monitoring the total arclength
uidx = coco_get_func_data(prob, 'dist', 'uidx');
prob = coco_add_pars(prob, 'length', uidx(end), 'L');
  
prob = coco_set(prob, 'cont', 'PtMX', [20 50]);

fprintf('\n Run=''%s'': Continue branch of isolas.\n', 'isola');

bd2  = coco(prob, 'isola', [], 1, {'L' 'lambda' 'tau'}, [0 1]);
  
%% Continue Hopf bifurcation curve

% The continuation problem encoded below includes two monitor functions
% that evaluate to each of the two problem parameters, respectively, and
% two corresponding inactive continuation parameters 'tau' and 'lambda'.
% Its dimensional deficit equals -1. The call to the coco entry-point
% function indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'lambda' and 'tau' are both released and allowed
% to vary during continuation.

prob = coco_prob();
prob = ode_HB2HB(prob, '', 'initial', labs(2));

fprintf(...
  '\n Run=''%s'': Continue HB points from point %d in run ''%s''.\n', ...
  'HB-curve', labs(2), 'initial');

bd3 = coco(prob, 'HB-curve', [], ... % run name = 'HB-curve'
  1, {'tau' 'lambda'}, [0 2]);       % continuation parameters and interval for 'tau'

%% Graphical representation

% Family of isolas and Hopf bifurcation curve
figure(1); clf; hold on; grid on; axis([0 1.5 0.2 1])
thm = struct();
thm.sol.RO = {'k.', 'MarkerSize', 5};
coco_plot_sol(thm, 'isola', 'isola', 1:N, 'p', 'x')
thm.special = {'BTP'};
coco_plot_bd(thm, 'HB-curve', 'tau', 'x')
hold off;

% Isola length against lambda
figure(2); clf
thm = struct('ylab', 'length');
coco_plot_bd(thm, 'isola', 'lambda', 'L', @(L) L*N);
grid on; axis([-inf inf -inf inf])
