function prob = ode_ep2HB(prob, oid, varargin)
%ODE_EP2HB   Start continuation of Hopf bifurcations of equilibrium points.
%
% PROB = ODE_EP2HB(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-ep-end' | '-end-ep' }
%
% Start or restart a continuation of Hopf bifurcations of equilibrium
% points from a previously obtained Hopf bifurcation. To start from a saved
% Hopf bifurcation, at least the name RUN of the continuation run and the
% solution label LAB must be given. The label LAB must be the label of a
% Hopf bifurcation.
%
% The arguments and their meaning are identical to ODE_EP2EP.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-ep-end' and '-end-ep' (optional). Either marks the end of input
%        to ODE_EP2HB.
% 
% See also: ODE_EP2EP, ODE_EP2SN, ODE_HB2HB, EP_READ_SOLUTION, EP_ADD,
% EP_ADD_HB

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_ep2HB.m 2897 2015-10-07 17:43:39Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-ep-end',       '',    '',    'end', {}
  '-end-ep',       '',    '',    'end', {}
  };
args = coco_parse(grammar, args_spec, opts_spec, varargin{:});

[sol, data] = ep_read_solution(args.soid, args.run, args.lab);
data = ode_init_data(prob, data, oid, 'ep');
[prob, data] = ep_add(prob, data, sol, '-no-test', '-cache-jac');
[prob, data] = ep_add_var(prob, data, sol.var.v);
prob = ep_add_HB(prob, data, sol);
prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info('HB'));

end
