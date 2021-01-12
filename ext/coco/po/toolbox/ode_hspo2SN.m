function prob = ode_hspo2SN(prob, oid, varargin)
%ODE_HSPO2SN   Start continuation of saddle-node bifurcations of multisegment periodic orbits.
%
% PROB = ODE_HSPO2SN(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-hspo-end' | '-end-hspo' }
%
% Start or restart a continuation of saddle-node bifurcations of
% multisegment periodic orbits from a previously obtained saddle-node
% bifurcation. To start from a saved saddle-node bifurcation, at least the
% name RUN of the continuation run and the solution label LAB must be
% given. The label LAB must be the label of a saddle-node bifurcation.
%
% The arguments and their meaning are identical to ODE_HSPO2HSPO.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-hspo-end' and '-end-hspo' (optional). Either marks the end of
%        input to ODE_HSPO2SN.
% 
% See also: ODE_HSPO2HSPO, ODE_HSPO2PD, ode_HSPO2TR, ODE_SN2SN,
% HSPO_READ_SOLUTION, HSPO_ADD, HSPO_ADD_SN

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_po2PD.m 2849 2015-05-17 20:32:46Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-hspo-end',       '',    '',    'end', {}
  '-end-hspo',       '',    '',    'end', {}
  };
args = coco_parse(grammar, args_spec, opts_spec, varargin{:});

sol  = hspo_read_solution(args.soid, args.run, args.lab);

stbid = coco_get_id(args.soid, 'hspo');
data  = coco_read_solution(stbid, args.run, args.lab, 'data');
data  = hspo_init_data(prob, data, oid, 'hspo_orb');

tsid = coco_get_id(oid, 'hspo.orb');
ssid = coco_get_id(stbid, 'orb');
prob = ode_bvp2bvp(prob, tsid, args.run, ssid, args.lab, ...
  '-var', sol.var.v);

[prob, data] = hspo_add(prob, data, '-no-test');
prob = hspo_add_SN(prob, data, sol);
prob = ode_add_tb_info(prob, oid, 'hspo', 'hspo', 'hspo', hspo_sol_info('SN'));

end
