function prob = adjt_isol2hspo(prob, oid, varargin)
%ADJT_ISOL2HSPO   Append adjoint of 'hspo' instance from initial guess.
%
% PROB = ADJT_ISOL2HSPO(PROB, OID, VARARGIN)
% VARARGIN = { [OPTS] }
% OPTS     = { '-hspo-end' | '-end-hspo'}
%
% Append adjoint of a 'hspo' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB using ODE_ISOL2HSPO. The preceding call to ODE_ISOL2HSPO must
% include explicit Jacobians, while functions evaluating second derivatives
% are optional.
%
% On input:
%
% PROB   : Continuation problem structure.
%
% OID    : Object instance identifier (string). The corresponding toolbox
%          instance identifier is coco_get_id(OID, 'po'). Pass the empty
%          string '' for a simple continuation of periodic orbits. Pass a
%          non-trivial object identifier if an instance of the HSPO toolbox
%          is part of a composite continuation problem.
%
% OPTS   : '-hspo-end' and '-end-hspo' (optional). Either '-hspo-end' or
%          '-end-hspo' marks the end of input to ADJT_ISOL2HSPO. 
%
% See also: ODE_ISOL2HSPO, HSPO_READ_ADJOINT, HSPO_ADJT_INIT_DATA,
% HSPO_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_isol2po.m 2928 2015-10-30 14:18:59Z hdankowicz $

tbid  = coco_get_id(oid, 'hspo');
fdata = coco_get_func_data(prob, tbid, 'data');

segoid = coco_get_id(tbid, 'orb');
prob   = adjt_isol2bvp(prob, segoid);

[sol, data] = hspo_read_adjoint('', '', fdata);
data = hspo_adjt_init_data(prob, data, oid, 'hspo_orb');
prob = hspo_construct_adjt(prob, data, sol);

end
