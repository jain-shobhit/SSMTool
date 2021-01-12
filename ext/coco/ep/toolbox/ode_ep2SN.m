function prob = ode_ep2SN(prob, oid, varargin)
%ODE_EP2SN   Start continuation of saddle-node bifurcations of equilibrium points.
%
% PROB = ODE_EP2SN(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-ep-end' | '-end-ep' }
%
% Start or restart a continuation of saddle-node bifurcations of
% equilibrium points from a previously obtained saddle-node bifurcation or
% fold point. To start from a saved sadle-node bifurcation or fold point,
% at least the name RUN of the continuation run and the solution label LAB
% must be given. The label LAB must be the label of a sadle-node
% bifurcation or fold point.
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
%        to ODE_EP2SN.
% 
% See also: ODE_EP2EP, ODE_EP2HB, ODE_SN2SN, EP_READ_SOLUTION, EP_ADD,
% EP_ADD_SN

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_ep2SN.m 2897 2015-10-07 17:43:39Z hdankowicz $

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
prob = ep_add_SN(prob, data, sol);
prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info('SN'));

end
