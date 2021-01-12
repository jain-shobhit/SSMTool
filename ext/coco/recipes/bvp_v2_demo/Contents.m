% Toolbox: 'bvp'
% Version:  2.0
% 
% Source:   Sect. 8.3 of Recipes for Continuation
%
% The example in bvp_v2_demo illustrates the value in allowing function
% data parameterizing boundary conditions to change during continuation.
%
% Scripts
%   <a href="matlab: coco_recipes_edit bvp_v2_demo demo_lienard">demo_lienard.m</a>  - Example 8.5.
%
% Functions
%   <a href="matlab: coco_recipes_edit bvp_v2_demo figure_8_7">figure_8_7.m</a>    - Generate figure 8.7.
%
% Other files
%   <a href="matlab: coco_recipes_edit bvp_v2_demo lienard">lienard.m</a>       - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit bvp_v2_demo per_bc">per_bc.m</a>        - 'bvp'-compatible encoding of boundary conditions.
%   <a href="matlab: coco_recipes_edit bvp_v2_demo per_bc_DFDX">per_bc_DFDX.m</a>   - 'bvp'-compatible encoding of Jacobian.
%   <a href="matlab: coco_recipes_edit bvp_v2_demo per_bc_update">per_bc_update.m</a> - Slot function to update instance data.
%
% See also <a href="matlab: coco_recipes_doc bvp_v2 bvp_v2">bvp_v2</a>, <a href="matlab: coco_recipes_doc coll_v1 coll_v1">coll_v1</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
