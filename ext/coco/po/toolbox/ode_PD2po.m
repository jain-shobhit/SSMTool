function prob = ode_PD2po(prob, oid, varargin)
%ODE_PD2PO   Start continuation of period-doubled periodic orbits from saved period-doubling bifurcation.
%
% PROB = ODE_PD2PO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-po-end' | '-end-po' | '-var' VECS }
%
% Start a continuation of periodic orbits of twice the period from a
% period-doubling bifurcation that was obtained and saved to disk in a
% previous continuation. To restart from a saved period-doubling
% bifurcation, at least the name RUN of the continuation run and the
% solution label LAB must be given.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string). Pass the empty string
%        '' for a simple continuation of periodic orbits.
%
% See ODE_ISOL2PO for more details on PROB and OID.
%
% RUN  : Run identifier (string or cell-array of strings). Name of the run
%        from which to restart a new continuation run.
% SOID : Source object instance identifier (string, optional). If the
%        argument SOID is omitted, OID is used. Pass the empty string ''
%        for OID and omit SOID for a simple continuation of multisegment
%        periodic orbits. Pass non-trivial object identifiers if an
%        instance of the PO toolbox is part of a composite continuation
%        problem.
% LAB  : Solution label (integer). The integer label assigned by COCO to an
%        trajectory segment during the continuation run RUN.
%
% OPTS : '-po-end', '-end-po', and '-var' VECS (optional, multiple options
%        may be given). Either '-po-end' or '-end-po' mark the end of input
%        to ODE_PD2PO. The option '-var' indicates the inclusion of the
%        variational problem for the corresponding trajectory segment,
%        where the initial solution guess for the perturbations to the
%        initial condition of the orbit is given by the content of VECS.
%
% See also: ODE_ISOL2PO, PO_READ_SOLUTION, ODE_PO2PO, PO_ADD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_PD2po.m 2947 2016-04-07 12:58:28Z hdankowicz $

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

[sol, data] = po_read_solution(args.soid, args.run, args.lab);
data = coco_read_solution(data.cid, args.run, args.lab, 'data');
cid = coco_get_id(oid,'po.orb.coll');
if isempty(coco_get(prob, cid, 'NTST'))
  prob = coco_set(prob, cid, 'NTST', 2*data.coll.NTST);
end
prob = ode_isol2po(prob, oid, data.fhan, data.dfdxhan, data.dfdphan, ...
  sol.pd_orb.t0, sol.pd_orb.x0, data.pnames, sol.pd_orb.p0, ...
  '-var', opts.vecs);

end
