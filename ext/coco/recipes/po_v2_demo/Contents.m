% Toolbox: 'po'
% Version:  2.0
% 
% Source:   Sect. 17.2 of Recipes for Continuation
%
% The examples in po_v2_demo illustrate the detection of bifurcations along
% a family of periodic orbits and branch switching at a period-doubling
% bifurcation point using chart data.
%
% Scripts
%   <a href="matlab: coco_recipes_edit po_v2_demo demo_po_v2">demo_po_v2.m</a> - Bifurcations and branch switching. [Not in Recipes for Continuation]
%
% Other files
%   <a href="matlab: coco_recipes_edit po_v2_demo tor">tor.m</a>        - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit po_v2_demo tor_DFDX">tor_DFDX.m</a>   - 'coll'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: coco_recipes_edit po_v2_demo tor_DFDP">tor_DFDP.m</a>   - 'coll'-compatible encoding of Jacobian w.r.t. problem parameters.
%
% See also <a href="matlab: coco_recipes_doc po_v2 po_v2">po_v2</a>, <a href="matlab: coco_recipes_doc varcoll_v1 varcoll_v1">varcoll_v1</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
