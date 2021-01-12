% Toolbox: 'atlas1d'
% Version:  7.0
% 
% Source:   Sect. 17.3.1 of Recipes for Continuation
%
% The example in atlas1d_v7_demo illustrates the use of automatic branch
% switching for bifurcation analysis of periodic orbits of a hybrid
% dynamical system.
%
% Scripts
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo demo_atlas1d_v7">demo_atlas1d_v7.m</a> - Examples 17.4 and 17.5.
%
% Functions
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo figure_17_6">figure_17_6.m</a>     - Generate figure 17.6.
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo figure_17_7">figure_17_7.m</a>     - Generate figure 17.7.
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo figure_17_8">figure_17_8.m</a>     - Generate figure 17.8.
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo figure_17_9">figure_17_9.m</a>     - Generate figure 17.9.
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo figure_17_10">figure_17_10.m</a>    - Generate figure 17.10.
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo figure_17_11">figure_17_11.m</a>    - Generate figure 17.11.
%
% Other files
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo duff">duff.m</a>            - 'hspo'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo duff_add_IP">duff_add_IP.m</a>     - Slot function for storing to bifurcation data.
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo duff_events">duff_events.m</a>     - 'hspo'-compatible encoding of event function.
%   <a href="matlab: coco_recipes_edit atlas1d_v7_demo duff_resets">duff_resets.m</a>     - 'hspo'-compatible encoding of reset function.
%
% See also <a href="matlab: coco_recipes_doc atlas1d_v7 atlas1d_v7">atlas1d_v7</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
