function prob = adjt_isol2bvp(prob, oid, varargin)
%ADJT_ISOL2BVP   Append adjoint of 'bvp' instance from initial guess.
%
% PROB     = ADJT_ISOL2BVP(PROB, OID, VARARGIN)
% VARARGIN = { [OPTS] }
% OPTS = { '-bvp-end' | '-end-bvp' }
%
% Append adjoint of a 'bvp' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB using ODE_ISOL2BVP. The preceding call to ODE_ISOL2BVP must include
% explicit Jacobians, while functions evaluating second derivatives are
% optional.
%
% On input:
%
% PROB       : Continuation problem structure.
%
% OID        : Object instance identifier (string). The corresponding
%              toolbox instance identifier is coco_get_id(OID, 'bvp'). Pass
%              the empty string '' for a simple continuation of collections
%              of constrained trajectory segments. Pass a non-trivial
%              object identifier if a 'bvp' instance is part of a composite
%              continuation problem.
%
% OPTS       : '-bvp-end' and '-end-bvp' (optional, multiple options may be
%               given). Either '-bvp-end' or '-end-bvp' marks the end of
%               input to ADJT_ISOL2BVP.
%
% See also: ODE_ISOL2BVP, BVP_READ_ADJOINT, BVP_ADJT_INIT_DATA,
% BVP_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_isol2bvp.m 2839 2015-03-05 17:09:01Z fschild $

tbid  = coco_get_id(oid, 'bvp');
fdata = coco_get_func_data(prob, tbid, 'data');

for i=1:fdata.nsegs
  segoid = coco_get_id(tbid, sprintf('seg%d', i));
  prob   = adjt_isol2coll(prob, segoid);
end

[sol, data] = bvp_read_adjoint('', '', fdata);
data = bvp_adjt_init_data(prob, data, oid, 'bvp_bc');
prob = bvp_construct_adjt(prob, data, sol);

end
