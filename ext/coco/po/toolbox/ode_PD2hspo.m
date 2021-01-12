function prob = ode_PD2hspo(prob, oid, varargin)
%ODE_PD2HSPO   Start continuation of period-doubled multisegment periodic orbits from saved period-doubling bifurcation.
%
% PROB = ODE_PD2HSPO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-hspo-end' | '-end-hspo' | '-var' VECS }
%
% Start a continuation of multisegment periodic orbits of twice the period
% from a period-doubling bifurcation that was obtained and saved to disk in
% a previous continuation. To restart from a saved period-doubling point,
% at least the name RUN of the continuation run and the solution label LAB
% must be given.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string). Pass the empty string
%        '' for a simple continuation of multisegment periodic orbits.
%
% See ODE_ISOL2HSPO for more details on PROB and OID.
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
% OPTS   : '-hspo-end', '-end-hspo', and '-var' VECS (optional, multiple
%          options may be given). Either '-hspo-end' or '-end-hspo' marks
%          the end of input to ODE_PDL2HSPO. The option '-var' indicates
%          the inclusion of the variational problem for each trajectory
%          segment, where the initial solution guess for the perturbations
%          to the initial conditions  of the i-th trajectory segment is
%          given by the content of the i-th element of the cell array VECS.
%
% See also: ODE_ISOL2HSPO, HSPO_READ_SOLUTION, ODE_HSPO2HSPO, HSPO_ADD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_PD2po.m 2863 2015-07-26 22:19:05Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-hspo-end',     '', '',  'end', {}
  '-end-hspo',     '', '',  'end', {}
       '-var', 'vecs', [], 'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

[sol, data] = hspo_read_solution(args.soid, args.run, args.lab);
nsegs = data.hspo_orb.nsegs;
for i=1:nsegs
  cid = coco_get_id(data.bvid, sprintf('seg%d.coll', i));
  cdata = coco_read_solution(cid, args.run, args.lab, 'data');
  cid = coco_get_id(oid, sprintf('hspo.orb.bvp.seg%d.coll', i));
  if isempty(coco_get(prob, cid, 'NTST'))
    prob = coco_set(prob, cid, 'NTST', cdata.coll.NTST);
  end
  cid = coco_get_id(oid, sprintf('hspo.orb.bvp.seg%d.coll', i+nsegs));
  if isempty(coco_get(prob, cid, 'NTST'))
    prob = coco_set(prob, cid, 'NTST', cdata.coll.NTST);
  end
end

prob = ode_isol2hspo(prob, oid, data.fhan, data.dfdxhan, data.dfdphan, ...
  sol.pd_orb.modes, sol.pd_orb.events, sol.pd_orb.resets, ...
  sol.pd_orb.t0, sol.pd_orb.x0, data.pnames, sol.pd_orb.p0, ...
  '-var', opts.vecs);

end
