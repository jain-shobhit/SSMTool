function prob = povar_isol2orb(prob, oid, varargin)
%POVAR_ISOL2ORB   Append 'po' instance and 'varcoll' instance constructed from initial data.
%
% Construct an instance of 'po' and append the corresponding variational
% zero problem.
%
% PROB     = POVAR_ISOL2ORB(PROB, OID, VARARGIN)
% VARARGIN = { PO @DFDXDX @DFDXDP}
%
% PROB    - Continuation problem structure.
% OID     - Object instance identifier (string).
% PO      - Argument to po_isol2orb.
% @DFDXDX - Function handle for vectorized encoding of 3d array of first
%           partials of vector-field Jacobian with respect to problem variables.
% @DFDXDP - Function handle for vectorized encoding of 3d array of first
%           partials of vector-field Jacobian with respect to problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: povar_isol2orb.m 2839 2015-03-05 17:09:01Z fschild $

str = coco_stream(varargin{:});
prob = po_isol2orb(prob, oid, str); % Create 'po' instance
oid  = coco_get_id(oid, 'po.seg');
dfdxdxhan = str.get;
dfdxdphan = str.get;
prob = var_coll_add(prob, oid, dfdxdxhan, dfdxdphan); % Append variational zero problem

end
