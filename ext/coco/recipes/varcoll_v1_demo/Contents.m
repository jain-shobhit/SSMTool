% Toolbox: 'varcoll'
% Version:  1.0
% 
% Source:   Sect. 10.1 of Recipes for Continuation
%
% The examples in varcoll_v1_demo illustrate the use of the <a href="matlab: coco_recipes_edit varcoll_v1 var_coll_add">var_coll_add</a>,
% <a href="matlab: coco_recipes_edit varcoll_v1 po_mult_add">po_mult_add</a>, and <a href="matlab: coco_recipes_edit varcoll_v1 hspo_mult_add">hspo_mult_add</a> toolbox constructors to include the
% computation of the fundamental solution of the variational problem for
% each 'coll' instance, as well as the Floquet multipliers for (hybrid)
% periodic orbits.
%
% Scripts
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo demo_linode">demo_linode.m</a>       - Example 10.2.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo demo_hybrid">demo_hybrid.m</a>       - Example 10.3.
%
% Other files
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo linode">linode.m</a>            - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo linode_DFDX">linode_DFDX.m</a>       - 'coll'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo linode_DFDP">linode_DFDP.m</a>       - 'coll'-compatible encoding of Jacobian w.r.t. problem parameters.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo lin_bc">lin_bc.m</a>            - 'bvp'-compatible encoding of boundary conditions.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo pwlin">pwlin.m</a>             - 'hspo'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo pwlin_DFDX">pwlin_DFDX.m</a>        - 'hspo'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo pwlin_DFDP">pwlin_DFDP.m</a>        - 'hspo'-compatible encoding of Jacobian w.r.t. problem parameters.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo pwlin_events">pwlin_events.m</a>      - 'hspo'-compatible encoding of event function.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo pwlin_resets">pwlin_resets.m</a>      - 'hspo'-compatible encoding of reset function.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo pwlin_events_DFDX">pwlin_events_DFDX.m</a> - 'hspo'-compatible encoding of Jacobian w.r.t states.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo pwlin_events_DFDP">pwlin_events_DFDP.m</a> - 'hspo'-compatible encoding of Jacobian w.r.t problem parameters.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo pwlin_resets_DFDX">pwlin_resets_DFDX.m</a> - 'hspo'-compatible encoding of Jacobian w.r.t states.
%   <a href="matlab: coco_recipes_edit varcoll_v1_demo pwlin_resets_DFDP">pwlin_resets_DFDP.m</a> - 'hspo'-compatible encoding of Jacobian w.r.t problem parameters.
%
% See also <a href="matlab: coco_recipes_doc varcoll_v1 varcoll_v1">varcoll_v1</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
