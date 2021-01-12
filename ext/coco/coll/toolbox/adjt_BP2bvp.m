function prob = adjt_BP2bvp(prob, oid, varargin)
%ADJT_BP2BVP   Append adjoint of 'bvp' instance at branch point.
%
% PROB     = ADJT_BP2BVP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-bvp-end' | '-end-bvp' }
%
% Append adjoint of a 'bvp' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB at the same branch point using ODE_BP2BVP. 
%
% The arguments and their meaning are identical to ADJT_BVP2BVP.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-bvp-end' and '-end-bvp' (optional, multiple options may be
%        given). Either '-bvp-end' or '-end-bvp' marks the
%        end of input to ADJT_BP2BVP.
%
% See also:ADJT_BVP2BVP, BVP_READ_ADJOINT, BVP_ADJT_INIT_DATA,
% BVP_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
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
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:}); %#ok<ASGLU>

[sol, data] = bvp_read_adjoint(args.soid, args.run, args.lab);
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
