% Toolbox: 'bvp'
% Version:  2.0
%
% Compatibility: 'coll', versions 1.0 and 2.0
% 
% Source:   Sect. 8.3 of Recipes for Continuation
%
% This version of the 'bvp' toolbox demonstrates the implementation of a
% boundary-value problem solver by applying boundary conditions to a single
% instance of the 'coll' toolbox representing a trajectory segment. The
% boundary conditions are parameterized by user data, which can be updated
% during continuation.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit bvp_v2 bvp_isol2seg">bvp_isol2seg</a>      - Append 'bvp' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit bvp_v2 bvp_sol2seg">bvp_sol2seg</a>       - Append 'bvp' instance constructed from saved data.
%   <a href="matlab:coco_recipes_edit bvp_v2 bvp_close_seg">bvp_close_seg</a>     - Append an instance of 'bvp' to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit bvp_v2 bvp_close_seg>bvp_F">bvp_close_seg>bvp_F</a>         - COCO-compatible wrapper to boundary
%                                 conditions.
%   <a href="matlab:coco_recipes_edit bvp_v2 bvp_close_seg>bvp_DFDU">bvp_close_seg>bvp_DFDU</a>      - COCO-compatible wrapper to linearization
%                                 of boundary conditions.
%   <a href="matlab:coco_recipes_edit bvp_v2 bvp_close_seg>bvp_bc_update">bvp_close_seg>bvp_bc_update</a> - COCO-compatible wrapper to boundary
%                                 condition function data update function.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit bvp_v2 bvp_arg_check">bvp_arg_check</a>     - Basic argument checking for 'bvp' toolbox. 
%                       [Change from Recipes for Continuation, 1st edition, page 186:
%                       omits unused prob argument.]
%   <a href="matlab:coco_recipes_edit bvp_v2 bvp_init_data">bvp_init_data</a>     - Initialize toolbox data for an instance of 'bvp'.
%   <a href="matlab:coco_recipes_edit bvp_v2 bvp_read_solution">bvp_read_solution</a> - Read 'bvp' solution and toolbox data from disk.
%
% Toolbox data (struct or func_data [contingent on bc_update])
%   fhan      : Function handle to encoding of boundary conditions.
%   dfdxhan   : Optional function handle to encoding of Jacobian w.r.t. T,
%               x0, x1, and p.
%   bc_data   : Optional boundary condition function data structure.
%   bc_update : Optional function handle to encoding for updating boundary
%               condition function data structure.
%   T_idx     : Index for interval length.
%   x0_idx    : Index array for trajectory end point at t=0.
%   x1_idx    : Index array for trajectory end point at t=1.
%   p_idx     : Index array for problem parameters.
%   tbid      : Toolbox instance identifier [contingent on bc_update].
%
% Solution (struct)
%   t       : Array of basepoint time instances.
%   x       : Array of trajectory basepoint values.
%   p       : Problem parameters.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit bvp_v2 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc bvp_v2_demo bvp_v2_demo">bvp_v2_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
