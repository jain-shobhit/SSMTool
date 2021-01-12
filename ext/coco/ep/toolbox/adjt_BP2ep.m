function prob = adjt_BP2ep(prob, oid, varargin)
%ADJT_BP2EP   Append adjoint of 'ep' instance at branch point.
%
% PROB = ADJT_BP2EP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-ep-end' | '-end-ep' }
%
% Append adjoint of an 'ep' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB at the same branch point using ODE_BP2EP. 
%
% The arguments and their meaning are identical to ADJT_EP2EP.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-ep-end' and '-end-ep' (optional, multiple options may be given).
%        Either '-ep-end' or '-end-ep' mark the end of input to ADJT_BP2EP.
% 
% See also: ADJT_EP2EP, EP_READ_ADJOINT, EP_ADJT_INIT_DATA,
% EP_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: adjt_BP2ep.m 2897 2015-10-07 17:43:39Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}  
  };
opts_spec = {
  '-ep-end',     '', '',  'end', {}
  '-end-ep',     '', '',  'end', {}
  };
args = coco_parse(grammar, args_spec, opts_spec, varargin{:}); 

[sol, data] = ep_read_adjoint(args.soid, args.run, args.lab);
data = ep_adjt_init_data(prob, data, oid);
prob = ep_construct_adjt(prob, data, sol);

end
