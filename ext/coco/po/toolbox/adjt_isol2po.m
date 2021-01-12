function prob = adjt_isol2po(prob, oid, varargin)
%ADJT_ISOL2PO   Append adjoint of 'po' instance from initial guess.
%
% PROB = ADJT_ISOL2PO(PROB, OID, VARARGIN)
% VARARGIN = { [OPTS] }
% OPTS     = { '-po-end' | '-end-po'}
%
% Append adjoint of a 'po' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB using ODE_ISOL2PO. The preceding call to ODE_ISOL2PO must include
% explicit Jacobians, while functions evaluating second derivatives are
% optional.
%
% On input:
%
% PROB   : Continuation problem structure.
%
% OID    : Object instance identifier (string). The corresponding toolbox
%          instance identifier is coco_get_id(OID, 'po'). Pass the empty
%          string '' for a simple continuation of periodic orbits. Pass a
%          non-trivial object identifier if an instance of the PO toolbox
%          is part of a composite continuation problem.
%
% OPTS   : '-po-end' and '-end-po' (optional). Either '-po-end' or
%          '-end-po' marks the end of input to ADJT_ISOL2PO. 
%
% See also: ODE_ISOL2PO, PO_READ_ADJOINT, PO_ADJT_INIT_DATA,
% PO_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_isol2po.m 2928 2015-10-30 14:18:59Z hdankowicz $

tbid  = coco_get_id(oid, 'po');
fdata = coco_get_func_data(prob, tbid, 'data');

segoid = coco_get_id(tbid, 'orb');
prob   = adjt_isol2coll(prob, segoid);

[sol, data] = po_read_adjoint('', '', fdata);
data = po_adjt_init_data(prob, data, oid, 'po_orb');
prob = po_construct_adjt(prob, data, sol);

end
