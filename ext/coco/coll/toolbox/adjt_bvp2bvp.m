function prob = adjt_bvp2bvp(prob, oid, varargin)
%ADJT_BP2BVP   Append adjoint of 'bvp' instance from saved solution.
%
% PROB     = ADJT_BP2BVP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-bvp-end' | '-end-bvp' }
%
% Append adjoint of a 'bvp' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB from the same saved solution using ODE_BVP2BVP.
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
%        for OID and omit SOID for a simple continuation of trajectory
%        segments. Pass non-trivial object identifiers if an instance of
%        the COLL toolbox is part of a composite continuation problem.
% LAB  : Solution label (integer). The integer label assigned by COCO to an
%        trajectory segment during the continuation run RUN.
%
% OPTS : '-switch', '-bvp-end', and '-end-bvp' (optional, multiple options
%        may be given). '-switch' instructs ADJT_BVP2BVP to switch branches
%        at the solution point, which must then be a branch point. Either
%        '-bvp-end' or '-end-bvp' marks the end of input to ADJT_BVP2BVP.
%
% See also: ADJT_ISOL2BVP, ADJT_BP2BVP, BVP_READ_ADJOINT,
% BVP_ADJT_INIT_DATA, BVP_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: msbvp_sol2segs.m 2839 2015-03-05 17:09:01Z fschild $

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
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:}); 

if opts.switch
  prob = adjt_BP2bvp(prob, oid, args.run, args.soid, args.lab);
  return
end

[sol, data] = bvp_read_adjoint(args.soid, args.run, args.lab);
if ~isempty(sol.tl0)
  % We restart at branch point but stay on original solution manifold =>
  % reset start direction to tangent vector of primary branch.
  sol.tl0 = sol.tl;
  sol.pars_tl0 = sol.pars_tl;
  for i=2:data.nsegs
    sfid = sprintf('shared%d', i-1);
    sol.([sfid, '_tl0']) = sol.([sfid, '_tl']);
  end
  if ~coco_exist('NullItMX', 'class_prop', prob, 'cont')
    % Compute improved tangent to new branch on same solution manifold to
    % allow change of parameters on restart at a branch-point.
    prob = coco_set(prob, 'cont', 'NullItMX', 1);
  end
end

data = bvp_adjt_init_data(prob, data, oid, 'bvp_bc');

stbid = coco_get_id(args.soid, 'bvp');
ttbid = coco_get_id(oid, 'bvp');
nsegs = data.nsegs;
cids  = cell(1,nsegs);
for i=1:nsegs
  toid = coco_get_id(ttbid, sprintf('seg%d', i));
  soid = coco_get_id(stbid, sprintf('seg%d', i));
  cids{i} = coco_get_id(toid, 'coll');
  prob = adjt_BP2coll(prob, toid, args.run, soid, args.lab);
end
data.cids = cids;

prob = bvp_construct_adjt(prob, data, sol);

end
