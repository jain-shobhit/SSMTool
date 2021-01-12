% Toolbox: 'hspo'
% Version:  1.0
% 
% Source:   Sects. 9.3.2 and 15.3.3 of Recipes for Continuation
%
% The examples in hspo_v1_demo illustrate the formulation of multisegment
% boundary-value problems for continuation of periodic orbits in
% autonomous, hybrid dynamical systems, including segments of different
% state-space dimension, as well as event detection associated with grazing
% contact with a discontinuity surface.
%
% Scripts
%   <a href="matlab: coco_recipes_edit hspo_v1_demo demo_pwlin">demo_pwlin.m</a>        - Example 9.2.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo demo_stickslip">demo_stickslip.m</a>    - Examples 9.3 and 9.4.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo demo_impact">demo_impact.m</a>       - Example 15.2.
%
% Functions
%   <a href="matlab: coco_recipes_edit hspo_v1_demo figure_9_3">figure_9_3.m</a>        - Generate figure 9.3.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo figure_9_4">figure_9_4.m</a>        - Generate figure 9.4.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo figure_9_5">figure_9_5.m</a>        - Generate figure 9.5.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo figure_15_4">figure_15_4.m</a>       - Generate figure 15.4.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo figure_15_5">figure_15_5.m</a>       - Generate figure 15.5.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo figure_15_6">figure_15_6.m</a>       - Generate figure 15.6.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo figure_15_7">figure_15_7.m</a>       - Generate figure 15.7.
%
% Other files
%   <a href="matlab: coco_recipes_edit hspo_v1_demo pwlin">pwlin.m</a>             - 'hspo'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo pwlin_DFDX">pwlin_DFDX.m</a>        - 'hspo'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo pwlin_DFDP">pwlin_DFDP.m</a>        - 'hspo'-compatible encoding of Jacobian w.r.t. problem parameters.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo pwlin_events">pwlin_events.m</a>      - 'hspo'-compatible encoding of event function.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo pwlin_events_DFDX">pwlin_events_DFDX.m</a> - 'hspo'-compatible encoding of Jacobian w.r.t states.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo pwlin_events_DFDP">pwlin_events_DFDP.m</a> - 'hspo'-compatible encoding of Jacobian w.r.t problem parameters.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo pwlin_resets">pwlin_resets.m</a>      - 'hspo'-compatible encoding of reset function.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo pwlin_resets_DFDX">pwlin_resets_DFDX.m</a> - 'hspo'-compatible encoding of Jacobian w.r.t states.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo pwlin_resets_DFDP">pwlin_resets_DFDP.m</a> - 'hspo'-compatible encoding of Jacobian w.r.t problem parameters.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo stickslip">stickslip.m</a>         - 'hspo'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo stickslip_DFDX">stickslip_DFDX.m</a>    - 'hspo'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo stickslip_events">stickslip_events.m</a>  - 'hspo'-compatible encoding of event function.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo stickslip_resets">stickslip_resets.m</a>  - 'hspo'-compatible encoding of reset function.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo impact">impact.m</a>            - 'hspo'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo impact_DFDX">impact_DFDX.m</a>       - 'hspo'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: coco_recipes_edit hspo_v1_demo impact_DFDP">impact_DFDP.m</a>       - 'hspo'-compatible encoding of Jacobian w.r.t. problem parameters.
%
% See also <a href="matlab: coco_recipes_doc hspo_v1 hspo_v1">hspo_v1</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
