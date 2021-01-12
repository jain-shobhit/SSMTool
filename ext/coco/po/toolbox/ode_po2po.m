function prob = ode_po2po(prob, oid, varargin)
%ODE_PO2PO   Start continuation of periodic orbits from saved solution.
%
% PROB = ODE_PO2PO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-po-end' | '-end-po' | '-var' VECS }
%
% Restart a continuation of periodic orbits from a periodic orbit that was
% obtained and saved to disk in a previous continuation. To restart from a
% saved periodic orbit, at least the name RUN of the continuation run and
% the solution label LAB must be given.
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
%        for OID and omit SOID for a simple continuation of periodic
%        orbits. Pass non-trivial object identifiers if an instance of the
%        PO toolbox is part of a composite continuation problem.
% LAB  : Solution label (integer). The integer label assigned by COCO to an
%        trajectory segment during the continuation run RUN.
%
% OPTS : '-switch', '-po-end', '-end-po', and '-var' VECS (optional,
%        multiple options may be given). '-switch' instructs ODE_PO2PO to
%        switch branches at the solution point, which must then be a branch
%        point. Either '-po-end' or '-end-po' mark the end of input to
%        ODE_PO2PO. The option '-var' indicates the inclusion of the
%        variational problem for the corresponding trajectory segment,
%        where the initial solution guess for the perturbations to the
%        initial condition of the orbit is given by the content of VECS.
%
% See also: ODE_ISOL2PO, PO_READ_SOLUTION, PO_ADD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_po2po.m 2928 2015-10-30 14:18:59Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
   'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
  'SOID',     '',   'str', 'soid', oid, 'read', {}
   'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-po-end',       '',    '',    'end', {}
  '-end-po',       '',    '',    'end', {}
  '-switch', 'switch', false, 'toggle', {}
     '-var',   'vecs',    [],   'read', {}
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
if opts.switch
  prob = ode_BP2coll(prob, tsid, args.run, ssid, args.lab, ...
    '-var', opts.vecs);
else
  prob = ode_coll2coll(prob, tsid, args.run, ssid, args.lab, ...
    '-var', opts.vecs);
end

if data.ode.autonomous
  prob = po_add(prob, data);
else
  prob = po_add(prob, data, '-no-phase');
end
prob = ode_add_tb_info(prob, oid, 'po', 'po', 'po', po_sol_info());

end
