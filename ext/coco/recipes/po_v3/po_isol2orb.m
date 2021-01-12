function prob = po_isol2orb(prob, oid, varargin)
%PO_ISOL2ORB   Append 'po' instance constructed from initial data.
%
% Construct an instance of 'coll' and append periodic boundary conditions
% and an integral phase condition.
%
% Identical to po_v1.
%
% PROB     = PO_ISOL2ORB(PROB, OID, VARARGIN)
% VARARGIN = { COLL }
%
% PROB   - Continuation problem structure.
% OID    - Object instance identifier (string).
% COLL   - Argument to coll_isol2seg.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_isol2orb.m 2839 2015-03-05 17:09:01Z fschild $

tbid   = coco_get_id(oid, 'po');   % Create toolbox instance identifier
str    = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
segoid = coco_get_id(tbid, 'seg'); % Create segment object instance identifier
prob   = coll_isol2seg(prob, segoid, str); % Construct 'coll' instance

data = struct();
data = po_init_data(prob, tbid, data); % Build toolbox data
prob = po_close_orb(prob, tbid, data); % Append boundary conditions

end
