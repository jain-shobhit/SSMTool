% Toolbox: 'atlas2d'
% Version:  6.0
% 
% Source:   Sects. 14.2 and 14.3 of Recipes for Continuation
%
% The examples in atlas2d_v6_demo illustrates the restriction of
% continuation to a computational domain defined by interval limits on
% active continuation parameters and support for starting continuation from
% a point on the boundary of the computational domain.
%
% Scripts
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo demo_cylinder">demo_cylinder.m</a> - Example 14.5.
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo demo_resonant">demo_resonant.m</a> - Figures 14.3 and 14.4a).
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo demo_arnold">demo_arnold.m</a>   - Figures 14.4b).
%
% Functions
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo figure_14_2">figure_14_2.m</a>   - Generate figure 14.2.
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo figure_14_3">figure_14_3.m</a>   - Generate figure 14.3.
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo figure_14_4">figure_14_4.m</a>   - Generate figure 14.4.
%
% Other files
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo cylinder">cylinder.m</a>      - COCO-compatible function encoding. [Overloads Matlab function]
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo lang">lang.m</a>          - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo lang_DFDX">lang_DFDX.m</a>     - 'coll'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo lang_DFDP">lang_DFDP.m</a>     - 'coll'-compatible encoding of Jacobian w.r.t. problem parameters.
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo po_bc">po_bc.m</a>         - 'bvp-compatible encoding of boundary conditions.
%   <a href="matlab: coco_recipes_edit atlas2d_v6_demo po_bc_DFDX">po_bc_DFDX.m</a>    - 'bvp-compatible encoding of Jacobian.
%
% See also <a href="matlab: coco_recipes_doc atlas2d_v6 atlas2d_v6">atlas2d_v6</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
