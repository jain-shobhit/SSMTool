% Toolbox: 'bvp'
% Version:  1.0
% 
% Source:   Sect. 8.1.1 of Recipes for Continuation
%
% The examples in bvp_v1_demo illustrate the use of the <a href="matlab: coco_recipes_edit bvp_v1 bvp_isol2seg">bvp_isol2seg</a> and
% <a href="matlab: coco_recipes_edit bvp_v1 bvp_sol2seg">bvp_sol2seg</a> toolbox constructors, the <a href="matlab: coco_recipes_edit bvp_v1 bvp_read_solution">bvp_read_solution</a> solution extractor,
% and the <a href="matlab: doc coco_bd_read">coco_bd_read</a>, <a href="matlab: doc coco_bd_col">coco_bd_col</a>, and <a href="matlab: doc coco_bd_labs">coco_bd_labs</a> core utilities.
%
% Scripts
%   <a href="matlab: coco_recipes_edit bvp_v1_demo demo_catenary">demo_catenary.m</a> - Example 8.1.
%   <a href="matlab: coco_recipes_edit bvp_v1_demo demo_bratu">demo_bratu.m</a>    - Example 8.2.
%
% Functions
%   <a href="matlab: coco_recipes_edit bvp_v1_demo figure_8_1">figure_8_1.m</a>    - Generate figure 8.1.
%   <a href="matlab: coco_recipes_edit bvp_v1_demo figure_8_2">figure_8_2.m</a>    - Generate figure 8.2.
%
% Other files
%   <a href="matlab: coco_recipes_edit bvp_v1_demo catn">catn.m</a>         - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit bvp_v1_demo catn_bc">catn_bc.m</a>      - 'bvp'-compatible encoding of boundary conditions.
%   <a href="matlab: coco_recipes_edit bvp_v1_demo catn_bc_DFDX">catn_bc_DFDX.m</a> - 'bvp'-compatible encoding of Jacobian.
%   <a href="matlab: coco_recipes_edit bvp_v1_demo brat">brat.m</a>         - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit bvp_v1_demo brat_bc">brat_bc.m</a>      - 'bvp'-compatible encoding of boundary conditions.
%   <a href="matlab: coco_recipes_edit bvp_v1_demo brat_bc_DFDX">brat_bc_DFDX.m</a> - 'bvp'-compatible encoding of Jacobian.
%
% See also <a href="matlab: coco_recipes_doc bvp_v1 bvp_v1">bvp_v1</a>, <a href="matlab: coco_recipes_doc coll_v1 coll_v1">coll_v1</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
