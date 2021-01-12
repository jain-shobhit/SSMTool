function [prob, data] = ep_add_var(prob, data, vecs)
%EP_ADD_VAR   Append 'ep' instance variational problem.
%
% [PROB DATA] = EP_ADD_VAR(PROB, DATA, VECS [OPTS])
% OPTS = {}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function appends the variational zero problem J*v=w with identifier
% coco_get_id(DATA.oid, 'ep.var') to the 'ep' instance with identifier
% coco_get_id(DATA.oid, 'ep'). The corresponding function data is stored in
% the field DATA.ep_var. The corresponding instance of the 'ep' equilibrium
% point zero problem must have been created with the option '-cache-jac'
% passed to EP_ADD. The variational zero problem uses the cached Jacobian
% J = DATA.ep_Jx.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'ep' toolbox instance data structure (see EP_ADD).
% VECS : Initial solution guess for array v.
%
% OPTS : No options implemented in the current version.
%
% On output:
%
% PROB : Continuation problem structure with 'ep.var' instance added.
% DATA : Toolbox data structure with added structure ep_var containing 
%        function data of the 'ep' variational zero problem. This structure
%        contains the following fields:
%
%        ep_var.nvec  : dimension of variational zero problem; number of
%            columns of v
%        ep_var.fid   : function identifier of 'ep' variational zero
%            problem; fid = coco_get_id(OID, 'ep.var');
%        ep_var.x_idx : context-independent indices of state variables in
%            vector u of continuation variables; x = u(ep_var.x_idx)
%        ep_var.p_idx : context-independent indices of problem parameters
%            in vector u of continuation variables; p = u(ep_var.p_idx)
%        ep_var.v_idx : context-independent indices of initial
%            perturbations in vector u of continuation variables;
%            v = u(ep_var.v_idx)
%        ep_var.w_idx : context-independent indices of final perturbations
%            in vector u of continuation variables; w = u(ep_var.w_idx)
%        ep_var.rows  : row index array used in sparse construction of
%            Jacobian of variational zero problem; not saved with
%            solution data
%        ep_var.cols  : column index array used in sparse construction of
%            Jacobian of variational zero problem; not saved with
%            solution data
%        ep_var.J     : preallocated space for Jacobian of variational
%            zero problem with respect to state variables and problem
%            parameters; not saved with solution data
%        ep_var.Jw    : Jacobian of variational zero problem with respect
%            to final perturbations; not saved with solution data
%
% See also: ODE_INIT_DATA, EP_ADD, EP_READ_SOLUTION, ODE_ISOL2EP,
% ODE_EP2EP, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_add_var.m 2899 2015-10-07 21:18:01Z hdankowicz $

[prob, data] = ep_var_construct_eqn(prob, data, vecs);
end
