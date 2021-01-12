function prob = ode_bvp2bvp(prob, oid, varargin)
%ODE_BVP2BVP   Start continuation of collections of constrained trajectory segments from saved solution.
%
% PROB     = ODE_BVP2BVP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-bvp-end' | '-end-bvp' | '-var' VECS }
%
% Restart a continuation of collections of constrained trajectory segments
% from a collection of constrained trajectory segments that was obtained
% and saved to disk in a previous continuation. To restart from a saved
% trajectory segment, at least the name RUN of the continuation run and the
% solution label LAB must be given.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string). Pass the empty string
%        '' for a simple continuation of trajectory segments.
%
% See ODE_ISOL2BVP for more details on PROB and OID.
%
% RUN  : Run identifier (string or cell-array of strings). Name of the run
%        from which to restart a new continuation run.
% SOID : Source object instance identifier (string, optional). If the
%        argument SOID is omitted, OID is used. Pass the empty string ''
%        for OID and omit SOID for a simple continuation of collections of
%        constrained trajectory segments. Pass non-trivial object
%        identifiers if an instance of the BVP toolbox is part of a
%        composite continuation problem.
% LAB  : Solution label (integer). The integer label assigned by COCO to an
%        trajectory segment during the continuation run RUN.
%
% OPTS : '-switch', '-bvp-end', '-end-bvp', and '-var' VECS (optional,
%        multiple options may be given). '-switch' instructs ODE_BVP2BVP to
%        switch branches at the solution point, which must then be a branch
%        point. Either '-bvp-end' or '-end-bvp' marks the end of input to
%        ODE_BVP2BVP. The option '-var' indicates the inclusion of a
%        variational problem for each trajectory segment, where the initial
%        solution guess for the perturbations to the initial conditions of
%        the i-th trajectory segment is given by the content of the i-th
%        element of the cell array VECS.
%
% See also: ODE_ISOL2BVP, BVP_READ_SOLUTION, ODE_COLL2COLL

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_bvp2bvp.m 2839 2015-03-05 17:09:01Z fschild $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-bvp-end',       '',    '',    'end', {}
  '-end-bvp',       '',    '',    'end', {}
   '-switch', 'switch', false, 'toggle', {}
      '-var',   'vecs',    [],   'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

stbid = coco_get_id(args.soid, 'bvp');
data  = coco_read_solution(stbid, args.run, args.lab);
data  = bvp_init_data(prob, data, oid);

ttbid = coco_get_id(oid, 'bvp');
nsegs = data.nsegs;
cids  = cell(1,nsegs);
if ~isempty(opts.vecs)
  assert(numel(opts.vecs)==nsegs, ...
    '%s: incompatible specification of arrays of perturbations', ...
    mfilename);
else
  opts.vecs = cell(1,nsegs);
end
for i=1:nsegs
  toid = coco_get_id(ttbid, sprintf('seg%d', i));
  soid = coco_get_id(stbid, sprintf('seg%d', i));
  cids{i} = coco_get_id(toid, 'coll');
  if opts.switch
    prob = ode_BP2coll(prob, toid, args.run, soid, args.lab, ...
      '-var', opts.vecs{i});
  else
    prob = ode_coll2coll(prob, toid, args.run, soid, args.lab, ...
      '-var', opts.vecs{i});
  end
end
data.cids = cids;

prob = bvp_close_segs(prob, data);
prob = ode_add_tb_info(prob, oid, 'bvp', 'segs', 'bvp', bvp_sol_info());

end
