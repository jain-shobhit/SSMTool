% Toolbox: 'coll'
% Version:  1.0
% 
% Source:   Sects. 7.3 and 18.1 of Recipes for Continuation
%
% The examples in coll_v1_demo illustrate the use of the <a href="matlab: coco_recipes_edit coll_v1 coll_isol2seg">coll_isol2seg</a>
% and <a href="matlab: coco_recipes_edit coll_v1 coll_sol2seg">coll_sol2seg</a> toolbox constructors, the <a href="matlab: coco_recipes_edit coll_v1 coll_read_solution">coll_read_solution</a> solution
% extractor, the <a href="matlab: doc coco_bd_read">coco_bd_read</a>, <a href="matlab: doc coco_bd_col">coco_bd_col</a>, and <a href="matlab: doc coco_bd_labs">coco_bd_labs</a> core
% utilities, and the dependence on the discretization order.
%
% Scripts
%   <a href="matlab: coco_recipes_edit coll_v1_demo demo_catenary">demo_catenary.m</a>    - Section 7.3.1.
%   <a href="matlab: coco_recipes_edit coll_v1_demo demo_huxley">demo_huxley.m</a>      - Section 7.3.2.
%   <a href="matlab: coco_recipes_edit coll_v1_demo demo_doedel">demo_doedel.m</a>      - Section 7.3.3.
%   <a href="matlab: coco_recipes_edit coll_v1_demo demo_pneta">demo_pneta.m</a>       - Example 18.1.
%
% Functions
%   <a href="matlab: coco_recipes_edit coll_v1_demo figure_7_1">figure_7_1.m</a>       - Generate figure 7.1.
%   <a href="matlab: coco_recipes_edit coll_v1_demo figure_7_2">figure_7_2.m</a>       - Generate figure 7.2.
%   <a href="matlab: coco_recipes_edit coll_v1_demo figure_7_3">figure_7_3.m</a>       - Generate figure 7.3.
%   <a href="matlab: coco_recipes_edit coll_v1_demo figure_7_4">figure_7_4.m</a>       - Generate figure 7.4.
%   <a href="matlab: coco_recipes_edit coll_v1_demo figure_18_1">figure_18_1.m</a>      - Generate figure 18.1.
%
% Other files
%   <a href="matlab: coco_recipes_edit coll_v1_demo catenary">catenary.m</a>         - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit coll_v1_demo huxley">huxley.m</a>           - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit coll_v1_demo huxley_bcs">huxley_bcs.m</a>       - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit coll_v1_demo huxley_isol2het">huxley_isol2het.m</a>  - Problem constructor from initial data.
%   <a href="matlab: coco_recipes_edit coll_v1_demo huxley_sol2het">huxley_sol2het.m</a>   - Problem constructor from stored data.
%   <a href="matlab: coco_recipes_edit coll_v1_demo huxley_close_het">huxley_close_het.m</a> - Glue together two 'coll' instances.
%   <a href="matlab: coco_recipes_edit coll_v1_demo doedel">doedel.m</a>           - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit coll_v1_demo doedel_bcs">doedel_bcs.m</a>       - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit coll_v1_demo doedel_evs">doedel_evs.m</a>       - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit coll_v1_demo doedel_isol2het">doedel_isol2het.m</a>  - Problem constructor from initial data.
%   <a href="matlab: coco_recipes_edit coll_v1_demo doedel_sol2het">doedel_sol2het.m</a>   - Problem constructor from stored data.
%   <a href="matlab: coco_recipes_edit coll_v1_demo doedel_close_het">doedel_close_het.m</a> - Glue together 'coll' and 'alg' instances.
%   <a href="matlab: coco_recipes_edit coll_v1_demo pneta">pneta.m</a>            - 'coll'-compatible encoding of vector field.
%
% See also <a href="matlab: coco_recipes_doc coll_v1 coll_v1">coll_v1</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
