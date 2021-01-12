% Toolbox: 'alg'
% Version:  7.0
% 
% Source:   Sect. 16.2.3 of Recipes for Continuation
%
% The examples in alg_v7_demo illustrate support for fold detection and
% location using embedded or nonembedded monitor functions.
%
% Scripts
%   <a href="matlab: coco_recipes_edit alg_v7_demo demo_alg_v7">demo_alg_v7.m</a> - The cusp normal form.
%
% Functions
%   <a href="matlab: coco_recipes_edit alg_v7_demo figure_16_4">figure_16_4.m</a> - Generate figure 16.4.
%   <a href="matlab: coco_recipes_edit alg_v7_demo figure_16_5">figure_16_5.m</a> - Generate figure 16.5.
%
% Other files
%   <a href="matlab: coco_recipes_edit alg_v7_demo cusp">cusp.m</a>        - 'alg'-compatible zero function.
%   <a href="matlab: coco_recipes_edit alg_v7_demo cusp_DFDX">cusp_DFDX.m</a>   - Jacobian w.r.t. problem variables.
%   <a href="matlab: coco_recipes_edit alg_v7_demo cusp_DFDP">cusp_DFDP.m</a>   - Jacobian w.r.t. problem parameters.
%
% See also <a href="matlab: coco_recipes_doc alg_v7 alg_v7">alg_v7</a>, <a href="matlab: coco_recipes_doc atlas2d_v6 atlas2d_v6">atlas2d_v6</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
