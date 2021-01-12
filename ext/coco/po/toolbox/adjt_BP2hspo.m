function prob = adjt_BP2hspo(prob, oid, varargin)
%ADJT_BP2HSPO   Append adjoint of 'hspo' instance at branch point.
%
% PROB     = ADJT_BP2HSPO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-hspo-end' | '-end-hspo' }
%
% Append adjoint of a 'hspo' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB at the same branch point using ODE_BP2HSPO. 
%
% The arguments and their meaning are identical to ADJT_HSPO2HSPO.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-hspo-end' and '-end-hspo' (optional, multiple options may be
%        given). Either '-hspo-end' or '-end-hspo' marks the end of input
%        to ADJT_BP2HSPO.
%
% See also:ADJT_HSPO2HSPO, HSPO_READ_ADJOINT, HSPO_ADJT_INIT_DATA,
% HSPO_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: msbvp_sol2segs.m 2839 2015-03-05 17:09:01Z fschild $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-hspo-end', '', '',  'end', {}
  '-end-hspo', '', '',  'end', {}
  };
args = coco_parse(grammar, args_spec, opts_spec, varargin{:});

[sol, data] = hspo_read_adjoint(args.soid, args.run, args.lab);
data = hspo_adjt_init_data(prob, data, oid, 'hspo_orb');

soid = coco_get_id(args.soid, 'hspo.orb');
toid = coco_get_id(oid, 'hspo.orb');
bvid = coco_get_id(toid, 'bvp');
prob = adjt_BP2bvp(prob, toid, args.run, soid, args.lab);
data.bvid = bvid;

prob = hspo_construct_adjt(prob, data, sol);

end
