function prob = adjt_BP2po(prob, oid, varargin)
%ADJT_BP2PO   Append adjoint of 'po' instance at branch point.
%
% PROB     = ADJT_BP2PO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-po-end' | '-end-po' }
%
% Append adjoint of a 'po' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB at the same branch point using ODE_BP2PO. 
%
% The arguments and their meaning are identical to ADJT_PO2PO.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-po-end' and '-end-po' (optional, multiple options may be
%        given). Either '-po-end' or '-end-po' marks the
%        end of input to ADJT_BP2PO.
%
% See also:ADJT_PO2PO, PO_READ_ADJOINT, PO_ADJT_INIT_DATA,
% PO_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: msbvp_sol2segs.m 2839 2015-03-05 17:09:01Z fschild $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-po-end', '', '',  'end', {}
  '-end-po', '', '',  'end', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:}); %#ok<ASGLU>

[sol, data] = po_read_adjoint(args.soid, args.run, args.lab);
data = po_adjt_init_data(prob, data, oid, 'po_orb');

soid = coco_get_id(args.soid, 'po.orb');
toid = coco_get_id(oid, 'po.orb');
cid  = coco_get_id(toid, 'coll');
prob = adjt_BP2coll(prob, toid, args.run, soid, args.lab);
data.cid = cid;

prob = po_construct_adjt(prob, data, sol);

end
