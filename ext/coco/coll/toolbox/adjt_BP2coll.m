function prob = adjt_BP2coll(prob, oid, varargin)
%ADJT_BP2COLL   Append adjoint of 'coll' instance at branch point.
%
% PROB = ADJT_BP2COLL(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-coll-end' | '-end-coll' }
%
% Append adjoint of a 'coll' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB at the same branch point using ODE_BP2COLL. 
%
% The arguments and their meaning are identical to ADJT_COLL2COLL.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-coll-end' and '-end-coll' (optional, multiple options may be
%        given). Either '-coll-end' or '-end-coll' mark the
%        end of input to ADJT_BP2COLL.
% 
% See also:ADJT_COLL2COLL, COLL_READ_ADJOINT, COLL_ADJT_INIT_DATA,
% COLL_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_BP2coll.m 2898 2015-10-07 21:17:13Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-coll-end',     '', '',  'end', {}
  '-end-coll',     '', '',  'end', {}
  };
args = coco_parse(grammar, args_spec, opts_spec, varargin{:}); 

tbid = coco_get_id(oid, 'coll');
data = coco_get_func_data(prob, tbid, 'data');

sol  = coll_read_adjoint(args.soid, args.run, args.lab);
data = coll_adjt_init_data(prob, data, oid, 'coll_seg');
prob = coll_construct_adjt(prob, data, sol);

end
