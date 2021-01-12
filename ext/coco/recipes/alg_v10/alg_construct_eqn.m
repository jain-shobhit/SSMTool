function prob = alg_construct_eqn(prob, tbid, data, sol)
%ALG_CONSTRUCT_EQN   Append an instance of 'alg' to problem.
%
% Add zero functions encoded in terms of problem variables and problem
% parameters, monitor functions that evaluate to the problem parameters,
% and corresponding inactive continuation parameters.
%
% Differs from alg_v9 by the use of an event handler to distinguish between
% Hopf bifurcations and neutral saddle points.
%
% PROB = ALG_CONSTRUCT_EQN(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_construct_eqn.m 2839 2015-03-05 17:09:01Z fschild $

prob = coco_add_func(prob, tbid, @alg_F, @alg_DFDU, data, 'zero', ...
  'u0', sol.u);
uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  prob = coco_add_pars(prob, fid, uidx(data.p_idx), data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');
if data.alg.norm % Optional slot function
  data.tbid = tbid;
  prob = coco_add_slot(prob, tbid, @alg_bddat, data, 'bddat');
end

switch data.alg.FO % Optional monitor function for fold detection
  case {'regular', 'active'}
    fid_FO = coco_get_id(tbid, 'test', 'FO');
    data.tbid = tbid;
    data = coco_func_data(data);
    prob = coco_add_func(prob, fid_FO, @alg_fold, ...
      @alg_fold_DFDU, data, data.alg.FO, fid_FO, ...
      'uidx', uidx, 'fdim', 1); % Function type given in data.alg.FO
    prob = coco_add_slot(prob, tbid, @alg_update, data, 'update'); % Update bordering vectors
    prob = coco_add_event(prob, 'FO', fid_FO, 0); % FO - event type
end

if data.alg.HB % Optional monitor function for Hopf bifurcation detection
  fid_HB = coco_get_id(tbid, 'test', 'HB');
  data.tfid = fid_HB;
  prob = coco_add_chart_data(prob, fid_HB, [], []); % Allocation of chart data
  prob = coco_add_func(prob, fid_HB, @alg_hopf, data, ...
    'regular', fid_HB, 'uidx', uidx, 'passChart'); % Chart is passed to monitor function
  prob = coco_add_event(prob, @alg_evhan_HB, data, 'SP', fid_HB, 0); % Event handler
end

end

function [data y] = alg_F(prob, data, u)
%ALG_F   COCO-compatible zero function wrapper.
%
% Supports zero functions encoded to expect two arrays of problem variables
% and problem parameters, respectively, as input arguments.
%
% Identical to alg_v1. 

x = u(data.x_idx); % Extract problem variables
p = u(data.p_idx); % Extract problem parameters

y = data.fhan(x, p);

end

function [data J] = alg_DFDU(prob, data, u)
%ALG_DFDU   COCO-compatible linearization of zero function wrapper.
%
% Supports linearization of zero functions encoded to expect two arrays of
% problem variables and problem parameters, respectively, as input
% arguments.
%
% Identical to alg_v7.

x  = u(data.x_idx); % Extract problem variables
p  = u(data.p_idx); % Extract problem parameters

J1 = alg_fhan_DFDX(data, x, p);
J2 = alg_fhan_DFDP(data, x, p);
J  = sparse([J1 J2]);

end

function [data y] = alg_fold(prob, data, u)
%ALG_FOLD   Fold monitor function.
%
% Return argument y equals 0 at point of vanishing determinant of Jacobian
% of vector field with respect to problem variables. Bordering vectors
% data.b and data.c are updated after each completed continuation step.
%
% Identical to alg_v7.
%
% see also <a href="matlab:coco_recipes_edit alg_v10 alg_construct_eqn alg_update">alg_construct_eqn>alg_update</a>.

x  = u(data.x_idx); % Extract problem variables
p  = u(data.p_idx); % Extract problem parameters
Jx = alg_fhan_DFDX(data, x, p);
v  = [Jx data.b; data.c' 0]\data.rhs; % data.rhs = [0 ... 0 1]'
y  = v(end);

end

function [data J] = alg_fold_DFDU(prob, data, u)
%ALG_FOLD_DFDU   Linearization of fold monitor function.
%
% Second-order Jacobians of the zero function are here approximated by
% a mid-point finite-difference scheme applied to the first-order
% directional derivatives on the right-hand sides below:
%
% phi_{xx}(x,p)[t,.]=d/dh phi_{x}(x+h*t,p) at h=0
% phi_{xp}(x,p)[t,.]=d/dh phi_{p}(x+h*t,p) at h=0
%
% Identical to alg_v7.

x  = u(data.x_idx); % Extract problem variables
p  = u(data.p_idx); % Extract problem parameters
Jx = alg_fhan_DFDX(data, x, p);
M  = [Jx data.b; data.c' 0];
v  = M\data.rhs;
w  = data.rhs'/M;

h  = 1.0e-4*(1+norm(x));
J0 = alg_fhan_DFDX(data, x-h*v(data.x_idx), p);
J1 = alg_fhan_DFDX(data, x+h*v(data.x_idx), p);
hx = -w(data.x_idx)*(0.5/h)*(J1-J0);

J0 = alg_fhan_DFDP(data, x-h*v(data.x_idx), p);
J1 = alg_fhan_DFDP(data, x+h*v(data.x_idx), p);
hp = -w(data.x_idx)*(0.5/h)*(J1-J0);

J  = [hx hp];

end

function data = alg_update(prob, data, cseg, varargin)
%ALG_UPDATE   Update bordering vectors of fold monitor function.
%
% Update algorithm seeks to ensure that data.b and data.c approximate the
% left and right singular vectors of unit length corresponding to the
% smallest singular value of the Jacobian of the zero function with respect
% to the problem variables.
%
% Identical to alg_v7.

chart  = cseg.src_chart; % Current chart
uidx   = coco_get_func_data(prob, data.tbid, 'uidx');
u      = chart.x(uidx);
x      = u(data.x_idx); % Extract problem variables
p      = u(data.p_idx); % Extract problem parameters
Jx     = alg_fhan_DFDX(data, x, p);
w      = data.rhs'/[Jx data.b; data.c' 0];
data.b = w(data.x_idx)';
data.b = data.b/norm(data.b);
v      = [Jx data.b; data.c' 0]\data.rhs;
data.c = v(data.x_idx);
data.c = data.c/norm(data.c);

end

function [data chart y] = alg_hopf(prob, data, chart, u)
%ALG_HOPF   Hopf monitor function.
%
% Return argument y equals 0 at point where the real parts of two
% eigenvalues of Jacobian of vector field with respect to problem variables
% add up to 0.
%
% Differs from alg_v9 by relying on chart data for eliminating redundant
% computation of the Jacobian w.r.t. the problem variables and the
% corresponding eigenvalues in the corresponding event handler.

cdata = coco_get_chart_data(chart, data.tfid); % Read chart data
if ~isempty(cdata) && isfield(cdata, 'la')
  la = cdata.la;
else
  x  = u(data.x_idx);
  p  = u(data.p_idx);
  Jx = alg_fhan_DFDX(data, x, p);
  la = eig(Jx);
  chart = coco_set_chart_data(chart, data.tfid, struct('la', la)); % Write chart data
end
la = la(data.la_idx1)+la(data.la_idx2); % Compute all sums of pairs of distinct eigenvalues
sc = abs(la);
y  = real(prod((2*la)./(max(1,sc)+sc)));

end

function [data cseg msg] = alg_evhan_HB(prob, data, cseg, cmd, msg)
%ALG_EVHAN_HB   Hopf bifurcation event handler.
%
% Distinguish between Hopf bifurcations (two complex conjugate, purely
% imaginary eigenvalues) and neutral saddle points (two equal in magnitude,
% opposite in sign, real eigenvalues).

switch cmd
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate') % After unsuccessful location
      msg.action = 'warn';
    else % After detection
      cdata = coco_get_chart_data(cseg.ptlist{1}, data.tfid);
      la0 = cdata.la;
      cdata = coco_get_chart_data(cseg.ptlist{end}, data.tfid);
      la1 = cdata.la;
      switch abs(sum(sign(real(la0)))-sum(sign(real(la1))))
        case 4 % Two eigenvalues cross the imaginary axis
          msg.point_type = 'HB';
          msg.action     = 'locate';
        case 0 % No eigenvalue crosses the imaginary axis
          msg.point_type = 'NSad';
          if data.alg.NSad % Optional detection
            msg.action   = 'locate';
          else
            msg.action   = 'finish';
          end
        otherwise
          msg.point_type = 'HB';
          msg.action     = 'warn';
          msg.wmsg       = 'could not determine type of event';
      end
      msg.idx = 1;
    end
  case 'check' % Add special point to curve segment
    msg.action = 'add';
    msg.finish = true;
end

end

function [data res] = alg_bddat(prob, data, command, varargin)
%ALG_BDDAT   Append to bifurcation data cell array.
%
% Function adds column containing Euclidean norm of array of
% problem variables to bifurcation data cell array.
%
% Identical to alg_v6.

res = {};
switch command
  case 'init'
    % The following line corrects a missing toolbox instance identifier in
    % Recipes for Continuation, 1st edition, page 108.
    res   = sprintf('||%s.x||', data.tbid); % Column header
  case 'data'
    chart = varargin{1}; % Current chart
    uidx  = coco_get_func_data(prob, data.tbid, 'uidx');
    u     = chart.x(uidx);
    res   = norm(u(data.x_idx),2); % Column data
end

end
