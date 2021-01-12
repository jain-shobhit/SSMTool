function prob = ode_po2PD(prob, oid, varargin)
%ODE_PO2PD   Start continuation of period-doubling bifurcations of periodic orbits.
%
% PROB = ODE_PO2PD(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-po-end' | '-end-po' }
%
% Start or restart a continuation of period-doubling bifurcations of
% periodic orbits from a previously obtained period-doubling bifurcation.
% To start from a saved period-doubling bifurcation, at least the name RUN
% of the continuation run and the solution label LAB must be given. The
% label LAB must be the label of a period-doubling bifurcation.
%
% The arguments and their meaning are identical to ODE_PO2PO.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-po-end' and '-end-po' (optional). Either marks the end of input
%        to ODE_PO2PD.
% 
% See also: ODE_PO2PO, ODE_PO2SN, ode_PO2TR, ODE_PD2PD, PO_READ_SOLUTION,
% PO_ADD, PO_ADD_PD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_po2PD.m 2928 2015-10-30 14:18:59Z hdankowicz $

grammar   = 'RUN [SOID] LAB OPTS';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-po-end',       '',    '',    'end', {}
  '-end-po',       '',    '',    'end', {}
  };
args = coco_parse(grammar, args_spec, opts_spec, varargin{:});

sol  = po_read_solution(args.soid, args.run, args.lab);

stbid = coco_get_id(args.soid, 'po');
data  = coco_read_solution(stbid, args.run, args.lab, 'data');
data  = po_init_data(prob, data, oid, 'ode');

tsid = coco_get_id(oid, 'po.orb');
ssid = coco_get_id(stbid, 'orb');
prob = ode_coll2coll(prob, tsid, args.run, ssid, args.lab, ...
  '-var', sol.var.v);

if data.ode.autonomous
  [prob, data] = po_add(prob, data, '-no-test');
else
  [prob, data] = po_add(prob, data, '-no-test', '-no-phase');
end
prob = po_add_PD(prob, data, sol);
prob = ode_add_tb_info(prob, oid, 'po', 'po', 'po', po_sol_info('PD'));

end
