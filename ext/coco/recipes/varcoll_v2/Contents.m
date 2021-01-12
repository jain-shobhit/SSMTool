% Toolbox: 'varcoll'
% Version:  2.0
% 
% Source:   Sect. 10.2 of Recipes for Continuation
%
% This toolbox demonstrates the implementation of a zero problem for
% approximating the solution of the first variational equation along a
% segment object.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit varcoll_v2 var_coll_add">var_coll_add</a>  - Append an instance of 'varcoll' to problem.
%   <a href="matlab:coco_recipes_edit varcoll_v2 po_mult_add">po_mult_add</a>   - Add slot function computing Floquet multipliers of
%                   single-segment periodic orbit.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit varcoll_v2 var_coll_add>var_coll_F">var_coll_add>var_coll_F</a> - Variational collocation zero problem.
%   <a href="matlab:coco_recipes_edit varcoll_v2 var_coll_add>var_coll_DFDU">var_coll_add>var_coll_DFDU</a> - Linearization of variational collocation zero problem.
%   <a href="matlab:coco_recipes_edit varcoll_v2 po_mult_add>po_mult_eigs_bddat">po_mult_add>po_mult_eigs_bddat</a>     - Add Floquet multipliers to
%                                        bifurcation data.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit varcoll_v2 var_coll_add>var_coll_init_data">var_coll_add>var_coll_init_data</a> - Initialize toolbox data for an
%                                     instance of 'varcoll'.
%   <a href="matlab:coco_recipes_edit varcoll_v2 var_coll_init_sol">var_coll_init_sol</a> - Build 'varcoll' initial solution guess.
%
% Toolbox data (struct)
%   coll_id   : 'coll' object instance identifier.
%   dfdxdxhan : Function handle for vectorized encoding of 3d array of
%               first partials of vector-field Jacobian with respect to
%               problem variables.
%   dfdxdphan : Function handle for vectorized encoding of 3d array of
%               first partials of vector-field Jacobian with respect to
%               problem parameters.
%   dim       : State-space dimension.
%   M1_idx    : Index array for Jacobian of time-T flow.
%   ubp_idx   : Index array for basepoint values of fundamental solution.
%   u_shp     : Shape for vectorization.
%   R         : Coefficient matrix.
%   Id        : Identity matrix.
%   jac       : Coefficient matrix.
%   dxdxrows1 : Row index array for vectorization.
%   dxdxrows2 : Row index array for vectorization.
%   dxdxrows3 : Row index array for vectorization.
%   dxdxcols1 : Column index array for vectorization.
%   dxdxcols2 : Column index array for vectorization.
%   dxdxcols3 : Column index array for vectorization.
%   dxdprows  : Row index array for vectorization.
%   dxdpcols  : Column index array for vectorization.
%   dxdp_shp  : Shape for vectorization.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit varcoll_v2 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc varcoll_v2_demo varcoll_v2_demo">varcoll_v2_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
