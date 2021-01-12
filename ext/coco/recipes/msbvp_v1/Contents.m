% Toolbox: 'msbvp'
% Version:  1.0
%
% Compatibility: 'coll', versions 1.0 and 2.0
% 
% Source:   Sect. 9.1 of Recipes for Continuation
%
% This version of the 'msbvp' toolbox demonstrates the implementation of a
% multi-segment boundary-value problem solver by applying boundary
% conditions to multiple instances of the 'coll' toolbox representing
% individual trajectory segments. The boundary conditions are parameterized
% by user data, which can be updated during continuation.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit msbvp_v1 msbvp_isol2seg">msbvp_isol2seg</a>      - Append 'msbvp' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit msbvp_v1 msbvp_sol2seg">msbvp_sol2seg</a>       - Append 'msbvp' instance constructed from saved data.
%   <a href="matlab:coco_recipes_edit msbvp_v1 msbvp_close_seg">msbvp_close_seg</a>     - Append an instance of 'msbvp' to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit msbvp_v1 msbvp_close_seg>msbvp_F">msbvp_close_seg>msbvp_F</a>         - COCO-compatible wrapper to boundary
%                                 conditions.
%   <a href="matlab:coco_recipes_edit msbvp_v1 msbvp_close_seg>msbvp_DFDU">msbvp_close_seg>msbvp_DFDU</a>      - COCO-compatible wrapper to linearization
%                                 of boundary conditions.
%   <a href="matlab:coco_recipes_edit msbvp_v1 msbvp_close_seg>msbvp_bc_update">msbvp_close_seg>msbvp_bc_update</a> - COCO-compatible wrapper to boundary
%                                 condition function data update function.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit msbvp_v1 msbvp_arg_check">msbvp_arg_check</a>     - Basic argument checking for 'msbvp' toolbox. 
%   <a href="matlab:coco_recipes_edit msbvp_v1 msbvp_init_data">msbvp_init_data</a>     - Initialize toolbox data for an instance of 'msbvp'.
%   <a href="matlab:coco_recipes_edit msbvp_v1 msbvp_read_solution">msbvp_read_solution</a> - Read 'msbvp' solution and toolbox data from disk.
%
% Toolbox data (struct or func_data [contingent on bc_update])
%   nsegs     : Number of 'coll' instances.
%   pnames    : Empty cell array or cell array of string labels for
%              continuation parameters tracking problem parameters.
%   fbchan    : Function handle to encoding of boundary conditions.
%   dfbcdxhan : Optional function handle to encoding of Jacobian w.r.t. T,
%               x0, x1, and p.
%   bc_data   : Optional boundary condition function data structure.
%   bc_update : Optional function handle to encoding for updating boundary
%               condition function data structure.
%   T_idx     : Index array for interval lengths.
%   x0_idx    : Index array for trajectory end points at t=0.
%   x1_idx    : Index array for trajectory end points at t=1.
%   p_idx     : Index array for problem parameters.
%   tbid      : Toolbox instance identifier [contingent on bc_update].
%
% Solution (array of structs)
%   t       : Array of basepoint time instances.
%   x       : Array of trajectory basepoint values.
%   p       : Problem parameters.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit msbvp_v1 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc msbvp_v1_demo msbvp_v1_demo">msbvp_v1_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
