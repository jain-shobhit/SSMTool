% Toolbox: 'msbvp'
% Version:  1.0
% 
% Source:   Sect. 9.2 of Recipes for Continuation
%
% The example in msbvp_v1_demo illustrates the use of all-to-all boundary
% conditions to couple a user-defined number of 'coll' instances
% representing trajectory segments on a quasiperiodic invariant torus.
%
% Scripts
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo demo_langford">demo_langford.m</a> - Example 9.1.
%
% Functions
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo figure_9_1">figure_9_1.m</a>    - Generate figure 9.1.
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo figure_9_2">figure_9_2.m</a>    - Generate figure 9.2.
%
% Other files
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo torus_bc">torus_bc.m</a>      - 'msbvp'-compatible encoding of boundary conditions.
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo torus_bc_DFDX">torus_bc_DFDX.m</a> - 'msbvp'-compatible encoding of Jacobian.
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo lang">lang.m</a>          - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo lang_DFDX">lang_DFDX.m</a>     - 'coll'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo lang_DFDP">lang_DFDP.m</a>     - 'coll'-compatible encoding of Jacobian w.r.t. problem parameters.
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo lang_red">lang_red.m</a>      - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit msbvp_v1_demo plot_tori">plot_tori.m</a>     - Visualization of computed invariant tori.
%
% See also <a href="matlab: coco_recipes_doc msbvp_v1 msbvp_v1">msbvp_v1</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
