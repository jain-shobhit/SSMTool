function prob = ode_BP2bvp(prob, oid, varargin)
%ODE_BP2BVP   Switch to secondary branch of collections of constrained trajectory segments at branch point.
%
% PROB     = ODE_BP2BVP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-bvp-end' | '-end-bvp' | '-var' VECS }
%
% Start a continuation of collections of constrained trajectory segments
% along a secondary branch intersecting a previously computed branch with
% name RUN in a branch point. To start from a saved branch point, at least
% the name RUN of the continuation run and the solution label LAB must be
% given. The label LAB must be the label of a branch point.
%
% The arguments and their meaning are identical to ODE_BVP2BVP.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : ''-bvp-end', '-end-bvp', and '-var' VECS (optional, multiple
%        options may be given). Either '-bvp-end' or '-end-bvp' marks the
%        end of input to ODE_BVP2BVP. The option '-var' indicates the
%        inclusion of a variational problem for each trajectory segment,
%        where the initial solution guess for the perturbations to the
%        initial conditions of the i-th trajectory segment is given by the
%        content of the i-th element of the cell array VECS.
%
% See also: ODE_BVP2BVP, BVP_READ_SOLUTION, ODE_BP2COLL

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: msbvp_sol2segs.m 2839 2015-03-05 17:09:01Z fschild $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-bvp-end',     '', '',  'end', {}
  '-end-bvp',     '', '',  'end', {}
      '-var', 'vecs', [], 'read', {}
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
  prob = ode_BP2coll(prob, toid, args.run, soid, args.lab, ...
     '-var', opts.vecs{i});
end
data.cids = cids;

prob = bvp_close_segs(prob, data);
prob = ode_add_tb_info(prob, oid, 'bvp', 'segs', 'bvp', bvp_sol_info());

end
