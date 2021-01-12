% Toolbox: 'alg'
% Version:  8.0
% 
% Source:   Sect. 16.3 of Recipes for Continuation
%
% The examples in alg_v8_demo illustrate support for fold continuation
% using the Moore-Spence formulation.
%
% Scripts
%   <a href="matlab: coco_recipes_edit alg_v8_demo demo_alg_v8">demo_alg_v8.m</a> - The cusp normal form. [Not in Recipes for Continuation]
%
% Other files
%   <a href="matlab: coco_recipes_edit alg_v8_demo cusp">cusp.m</a>        - 'alg'-compatible zero function.
%   <a href="matlab: coco_recipes_edit alg_v8_demo cusp_DFDX">cusp_DFDX.m</a>   - Jacobian w.r.t. problem variables.
%   <a href="matlab: coco_recipes_edit alg_v8_demo cusp_DFDP">cusp_DFDP.m</a>   - Jacobian w.r.t. problem parameters.
%
% See also <a href="matlab: coco_recipes_doc alg_v8 alg_v8">alg_v8</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
