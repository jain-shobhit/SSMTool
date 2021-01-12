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
%   <a href="matlab: edit demo_atlas2d_v6">demo_atlas2d_v6.m</a>
%   <a href="matlab: edit demo_resonant">demo_resonant.m</a> - Figures 14.3 and 14.4a).
%
% Other files
%   <a href="matlab: edit torus">torus.m</a>          - 'coco'-compatible encoding of zero function.
%   <a href="matlab: edit angles">angles.m</a>          - 'coco'-compatible encoding of zero function.
%   <a href="matlab: edit lang">lang.m</a>          - 'coll'-compatible encoding of vector field.
%   <a href="matlab: edit lang_DFDX">lang_DFDX.m</a>     - 'coll'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: edit lang_DFDP">lang_DFDP.m</a>     - 'coll'-compatible encoding of Jacobian w.r.t. problem parameters.
%   <a href="matlab: edit po_bc">po_bc.m</a>         - 'bvp-compatible encoding of boundary conditions.
%   <a href="matlab: edit po_bc_DFDX">po_bc_DFDX.m</a>    - 'bvp-compatible encoding of Jacobian.
%
% See also <a href="matlab: doc ../atlas2d_v6/Contents.m">atlas2d_v6</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 3126 2019-06-12 02:46:20Z hdankowicz $
