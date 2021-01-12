% Toolbox: 'bvp'
% Version:  1.0
%
% Compatibility: 'coll', versions 1.0 and 2.0
% 
% Source:   Sect. 8.1.1 of Recipes for Continuation
%
% This version of the 'bvp' toolbox demonstrates the implementation of a
% boundary-value problem solver by applying boundary conditions to a single
% instance of the 'coll' toolbox representing a trajectory segment.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit bvp_v1 bvp_isol2seg">bvp_isol2seg</a>      - Append 'bvp' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit bvp_v1 bvp_sol2seg">bvp_sol2seg</a>       - Append 'bvp' instance constructed from saved data.
%   <a href="matlab:coco_recipes_edit bvp_v1 bvp_close_seg">bvp_close_seg</a>     - Append an instance of 'bvp' to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit bvp_v1 bvp_close_seg>bvp_F">bvp_close_seg>bvp_F</a>    - COCO-compatible wrapper to boundary conditions.
%   <a href="matlab:coco_recipes_edit bvp_v1 bvp_close_seg>bvp_DFDU">bvp_close_seg>bvp_DFDU</a> - COCO-compatible wrapper to linearization of
%                            boundary conditions.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit bvp_v1 bvp_arg_check">bvp_arg_check</a>     - Basic argument checking for 'bvp' toolbox. 
%                       [Change from Recipes for Continuation, 1st edition, page 186:
%                       omits unused prob argument.]
%   <a href="matlab:coco_recipes_edit bvp_v1 bvp_init_data">bvp_init_data</a>     - Initialize toolbox data for an instance of 'bvp'.
%   <a href="matlab:coco_recipes_edit bvp_v1 bvp_read_solution">bvp_read_solution</a> - Read 'bvp' solution and toolbox data from disk.
%
% Toolbox data (struct)
%   fhan    : Function handle to encoding of boundary conditions.
%   dfdxhan : Optional function handle to encoding of Jacobian w.r.t. T,
%             x0, x1, and p.
%   T_idx   : Index for interval length.
%   x0_idx  : Index array for trajectory end point at t=0.
%   x1_idx  : Index array for trajectory end point at t=1.
%   p_idx   : Index array for problem parameters.
%
% Solution (struct)
%   t       : Array of basepoint time instances.
%   x       : Array of trajectory basepoint values.
%   p       : Problem parameters.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit bvp_v1 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc bvp_v1_demo bvp_v1_demo">bvp_v1_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
