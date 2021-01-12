% Toolbox: 'hspo'
% Version:  2.0
% 
% Source:   Sect. 17.2.2 and 17.3.3 of Recipes for Continuation
%
% The examples in hspo_v2_demo illustrate the bifurcation analysis of
% periodic orbits of hybrid dynamical systems, including branch switching
% at a period-doubling point.
%
% Scripts
%   <a href="matlab: coco_recipes_edit hspo_v2_demo demo_hspo_v2">demo_hspo_v2.m</a> - Examples 17.3 and 17.5.
%
% Functions
%   <a href="matlab: coco_recipes_edit hspo_v2_demo figure_17_4">figure_17_4.m</a>  - Generate figure 17.4.
%   <a href="matlab: coco_recipes_edit hspo_v2_demo figure_17_5">figure_17_5.m</a>  - Generate figure 17.5.
%
% Other files
%   <a href="matlab: coco_recipes_edit hspo_v2_demo duff">duff.m</a>         - 'hspo'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit hspo_v2_demo duff_add_IP">duff_add_IP.m</a>  - Slot function for storing to bifurcation data.
%   <a href="matlab: coco_recipes_edit hspo_v2_demo duff_events">duff_events.m</a>  - 'hspo'-compatible encoding of event function.
%   <a href="matlab: coco_recipes_edit hspo_v2_demo duff_resets">duff_resets.m</a>  - 'hspo'-compatible encoding of reset function.
%
% See also <a href="matlab: coco_recipes_doc hspo_v2 hspo_v2">hspo_v2</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
