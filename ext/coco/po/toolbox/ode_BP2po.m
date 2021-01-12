function prob = ode_BP2po(prob, oid, varargin)
%ODE_BP2PO   Switch to secondary branch of periodic orbits at branch point.
%
% PROB = ODE_BP2PO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-po-end' | '-end-po' | '-var' VECS }
%
% Start a continuation of periodic orbits along a secondary branch
% intersecting a previously computed branch with name RUN in a branch
% point. To start from a saved branch point, at least the name RUN of the
% continuation run and the solution label LAB must be given. The label LAB
% must be the label of a branch point.
%
% The arguments and their meaning are identical to ODE_PO2PO.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-po-end', '-end-po', and '-var' VECS (optional, multiple options
%        may be given). Either '-po-end' or '-end-po' mark the end of input
%        to ODE_BP2PO. The option '-var' indicates the inclusion of the
%        variational problem for the corresponding trajectory segment,
%        where the initial solution guess for the perturbations to the
%        initial condition of the orbit is given by the content of VECS.
% 
% See also: ODE_PO2PO, PO_READ_SOLUTION, PO_ADD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_BP2po.m 2928 2015-10-30 14:18:59Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-po-end',     '', '',  'end', {}
  '-end-po',     '', '',  'end', {}
     '-var', 'vecs', [], 'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

stbid = coco_get_id(args.soid, 'po');
data  = coco_read_solution(stbid, args.run, args.lab, 'data');
data  = po_init_data(prob, data, oid, 'ode');
if data.po.bifus
    prob = coco_set(prob, data.cid, 'var', true);
end
tsid = coco_get_id(oid, 'po.orb');
ssid = coco_get_id(stbid, 'orb');
prob = ode_BP2coll(prob, tsid, args.run, ssid, args.lab, ...
  '-var', opts.vecs);

if data.ode.autonomous
  prob = po_add(prob, data);
else
  prob = po_add(prob, data, '-no-phase');
end
prob = ode_add_tb_info(prob, oid, 'po', 'po', 'po', po_sol_info());

end
