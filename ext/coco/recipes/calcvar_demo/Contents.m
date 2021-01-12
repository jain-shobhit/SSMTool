% Calculus of variations
% 
% Source:   Sect. 1.3 of Recipes for Continuation
%
% The examples in calcvar_demo illustrates the formulation of a system of
% algebraic equations, a discretized two-point boundary value problem, and
% a discretized optimization problem for continuation a family of extremal
% functions of the catenary problem from the calculus of variations.
%
% Scripts
%   <a href="matlab: coco_recipes_edit calcvar_demo demo_algebraic">demo_algebraic.m</a>    - Figure 1.5.
%   <a href="matlab: coco_recipes_edit calcvar_demo demo_differential">demo_differential.m</a> - Figure 1.8a).
%   <a href="matlab: coco_recipes_edit calcvar_demo demo_quadrature">demo_quadrature.m</a>   - Figure 1.8b).
%
% Functions
%   <a href="matlab: coco_recipes_edit calcvar_demo figure_1_5">figure_1_5.m</a>        - Generate figure 1.5.
%   <a href="matlab: coco_recipes_edit calcvar_demo figure_1_6">figure_1_6.m</a>        - Generate figure 1.6.
%   <a href="matlab: coco_recipes_edit calcvar_demo figure_1_7">figure_1_7.m</a>        - Generate figure 1.7.
%   <a href="matlab: coco_recipes_edit calcvar_demo figure_1_8">figure_1_8.m</a>        - Generate figure 1.8.
%
% Other files
%   <a href="matlab: coco_recipes_edit calcvar_demo caty_alg">caty_alg.m</a>          - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit calcvar_demo caty_ode">caty_ode.m</a>          - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit calcvar_demo caty_calcvar">caty_calcvar.m</a>      - 'calcvar'-compatible encoding of derivatives of Lagrangian.
%   <a href="matlab: coco_recipes_edit calcvar_demo calcvar_F">calcvar_F.m</a>         - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit calcvar_demo calcvar_DFDU">calcvar_DFDU.m</a>      - COCO-compatible encoding of Jacobian.
%   <a href="matlab: coco_recipes_edit calcvar_demo calcvar_system">calcvar_system.m</a>    - Utility for initialization of 'calcvar' function data.
%
% See also <a href="matlab: coco_recipes_doc coll_v1 coll_v1">coll_v1</a>, <a href="matlab: coco_recipes_doc bvp_v1 bvp_v1">bvp_v1</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
