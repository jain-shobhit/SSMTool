function prob = adjt_isol2coll(prob, oid, varargin)
%ADJT_ISOL2COLL   Append adjoint of 'coll' instance from initial guess.
%
% PROB = ADJT_ISOL2COLL(PROB, OID, VARARGIN)
% VARARGIN = { [OPTS] }
% OPTS = { '-coll-end' | '-end-coll' }
%
% Append adjoint of a 'coll' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB using ODE_ISOL2COLL. The preceding call to ODE_ISOL2COLL must
% include explicit Jacobians, while functions evaluating second derivatives
% are optional. 
%
% On input:
%
% PROB   : Continuation problem structure.
%
% OID    : Object instance identifier (string). The corresponding toolbox
%          instance identifier is coco_get_id(OID, 'coll'). Pass the empty
%          string '' for a simple continuation of trajectory segments. Pass
%          a non-trivial object identifier if an instance of the COLL
%          toolbox is part of a composite continuation problem.
%
% OPTS   : '-coll-end' and '-end-coll' (optional, multiple options may be
%          given). Either '-coll-end' or '-end-coll' marks the end of input
%          to ADJT_ISOL2COLL.
%
% See also: ODE_ISOL2COLL, COLL_READ_ADJOINT, COLL_ADJT_INIT_DATA,
% COLL_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_isol2coll.m 2898 2015-10-07 21:17:13Z hdankowicz $

tbid = coco_get_id(oid, 'coll');
data = coco_get_func_data(prob, tbid, 'data');

sol  = coll_read_adjoint('', '', data);
data = coll_adjt_init_data(prob, data, oid, 'coll_seg');
prob = coll_construct_adjt(prob, data, sol);

end
