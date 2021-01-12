function [prob, data] = coll_add_var(prob, data, vecs)
%COLL_ADD_VAR   Append 'coll' instance variational problem.
%
% [PROB DATA] = COLL_ADD_VAR(PROB, DATA, VECS [OPTS])
% OPTS = {}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function appends the variational problem J*v=0 with identifier
% coco_get_id(DATA.oid, 'coll.var') to the 'coll' instance with identifier
% coco_get_id(DATA.oid, 'coll'). The corresponding function data is stored
% in the field DATA.coll_var. The corresponding instance of the 'coll'
% trajectory segment zero problem at node coco_get_id(DATA.oid, 'coll')
% must have been created with the option '-cache-jac' passed to COLL_ADD.
% The variational zero problem uses the cached vector-field Jacobian J =
% DATA.coll_Jx.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'coll' toolbox instance data structure (see COLL_ADD).
% VECS : Initial solution guess for array of perturbations to initial point on trajectory.
%
% OPTS : No options implemented in the current version.
%
% On output:
%
% PROB : Continuation problem structure with with 'coll.var' instance added.
% DATA : Toolbox data structure with added structure coll_var containing 
%        function data of the 'coll' variational zero problem. This
%        structure contains the following fields:
%
%        coll_var.nvec    : dimension of variational zero problem; number of
%            columns of v
%        coll_var.fid     : function identifier of 'coll' variational zero
%            problem; fid = coco_get_id(OID, 'coll.var');
%        coll_var.xbp_idx : context-independent indices of base point
%            values of trajectory segment in vector u of continuation
%            variables; xbp = u(coll_var.xbp_idx)
%        coll_var.T_idx   : context-independent indices of interval
%            duration in vector u of continuation variables; T =
%            u(coll_var.T_idx)
%        coll_var.p_idx   : context-independent indices of problem
%            parameters in vector u of continuation variables; p =
%            u(coll_var.p_idx)
%        coll_var.v_idx   : context-independent indices of base point
%            values of solution to variational problem in vector u of
%            continuation variables; v = u(coll_var.v_idx)
%        coll_var.v0_idx  : context-independent indices of perturbation to
%            initial trajecory end point in vector u of continuation
%            variables; v0 = u(coll_var.v0_idx)
%        coll_var.v1_idx  : context-independent indices of perturbation to
%            final trajecory end point in vector u of continuation
%            variables; v1 = u(coll_var.v1_idx)
%        coll_var.tbp and coll_var.tbp_idx : temporal mesh and index array
%            used in remeshing variational problem; not saved with
%            solution data
%        coll_var.xtr, coll_var.xtrbeg, and coll_var.xtrend : index arrays
%            used in remeshing variational problem; not saved with
%            solution data
%        coll_var.xbp_shp, coll_var.rows, coll_var.cols, coll_var.J,
%        coll_var.Q, coll_var.Qv, and coll_var.zeros : shape and index
%            arrays, and preallocated matrix content for constructing
%            variational zero problem and its Jacobian; not saved with
%            solution data
%
% See also: ODE_INIT_DATA, COLL_ADD, COLL_READ_SOLUTION, ODE_ISOL2COLL,
% ODE_COLL2COLL, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_add_var.m 2898 2015-10-07 21:17:13Z hdankowicz $

[prob, data] = coll_var_construct_eqn(prob, data, vecs);
end
