% Toolbox: 'po'
% Version:  3.0
% 
% Source:   Sect. 20.3.1 of Recipes for Continuation
%
% The examples in po_v3_demo illustrate the use of adaptive remeshing using
% fixed- and variable-order moving mesh algorithms in the continuation of a
% family of approximate homoclinic trajectories.
% 
% Scripts
%   <a href="matlab: coco_recipes_edit po_v3_demo demo_po_v3">demo_po_v3.m</a>   - Figures 20.7-13.
%
% Functions
%   <a href="matlab: coco_recipes_edit po_v3_demo figure_20_7">figure_20_7.m</a>  - Generate figure 20.7.
%   <a href="matlab: coco_recipes_edit po_v3_demo figure_20_8">figure_20_8.m</a>  - Generate figure 20.8.
%   <a href="matlab: coco_recipes_edit po_v3_demo figure_20_9">figure_20_9.m</a>  - Generate figure 20.9.
%   <a href="matlab: coco_recipes_edit po_v3_demo figure_20_10">figure_20_10.m</a> - Generate figure 20.10.
%   <a href="matlab: coco_recipes_edit po_v3_demo figure_20_11">figure_20_11.m</a> - Generate figure 20.11.
%   <a href="matlab: coco_recipes_edit po_v3_demo figure_20_12">figure_20_12.m</a> - Generate figure 20.12.
%   <a href="matlab: coco_recipes_edit po_v3_demo figure_20_13">figure_20_13.m</a> - Generate figure 20.13.
%
% Other files
%   <a href="matlab: coco_recipes_edit po_v3_demo marsden">marsden.m</a>      - 'coll'-compatible encoding of vector field.
%
% See also <a href="matlab: coco_recipes_doc po_v3 po_v3">po_v3</a>, <a href="matlab: coco_recipes_doc po_v1_demo po_v1_demo">po_v1_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
