% Toolbox: 'coll'
% Version:  2.0
% 
% Source:   Sect. 9.3.2 of Recipes for Continuation
%
% The example in coll_v2_demo illustrates the construction of a segment
% continuation problem object using a zero-length initial solution guess.
%
% Scripts
%   <a href="matlab: coco_recipes_edit coll_v2_demo demo_coll_v2">demo_coll_v2.m</a> - Partial end-point conditions. [Not in Recipes for Continuation]
%
% Other files
%   <a href="matlab: coco_recipes_edit coll_v2_demo linode">linode.m</a>       - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit coll_v2_demo linode_DFDX">linode_DFDX.m</a>  - 'coll'-compatible encoding of Jacobian w.r.t. state.
%   <a href="matlab: coco_recipes_edit coll_v2_demo linode_DFDP">linode_DFDP.m</a>  - 'coll'-compatible encoding of Jacobian w.r.t. parameters.
%   <a href="matlab: coco_recipes_edit coll_v2_demo lin_bc">lin_bc.m</a>       - 'bvp'-compatible encoding of boundary conditions.
%   <a href="matlab: coco_recipes_edit coll_v2_demo add_IP">add_IP.m</a>       - Slot function for storing to bifurcation data.
%
% See also <a href="matlab: coco_recipes_doc coll_v2 coll_v2">coll_v2</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
