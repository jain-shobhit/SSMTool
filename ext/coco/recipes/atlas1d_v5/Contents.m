% Toolbox: 'atlas1d'
% Version:  5.0
% 
% Source:   Sect. 14.1 of Recipes for Continuation
%
% This toolbox implements a 1-dimensional expanding boundary atlas
% algorithm with constant step size and theta-based projection condition
% and predictor. The algorithm may start from a computational boundary and
% performs continuation in each direction along the solution manifold until
% the boundary of the computational domain is reached or the maximum number
% of steps is reached.
%
% Class definition
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min atlas_1d_min">atlas_1d_min</a> - 1d atlas class.
%
% Interface methods
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min init_prcond">init_prcond</a>     - Initialize projection condition.
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min init_chart">init_chart</a>      - Initialize chart.
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min init_admissible">init_admissible</a> - Remove inadmissible directions.
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min init_atlas">init_atlas</a>      - Initialize atlas.
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min flush">flush</a>           - Flush point list.
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min predict">predict</a>         - Compute predictor.
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min add_chart">add_chart</a>       - Add chart to point list.
%
% Private properties
%   boundary     : Cell array to hold atlas boundary charts.
%   cont         : Struct used to hold class settings.
%
% Private methods
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min get_settings">get_settings</a> - Defines default class settings.
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min merge">merge</a> - Merge curve segment into atlas.
%   <a href="matlab: coco_recipes_edit atlas1d_v5/@atlas_1d_min isneighbor">isneighbor</a> - Check if two charts are neighbors.
%
% Public methods
%   <a href="matlab:coco_recipes_edit atlas1d_v5/@atlas_1d_min atlas_1d_min>create">create</a> - Static construction method.
%
% Class settings
%   h     : Step size (default: 0.1)
%   PtMX  : Maximum number of continuation steps (default: 50)
%   theta : Theta method (default: 0.5)
%   almax : Critical angle between successive tangent vectors (default: 10)
%   Rmarg : Margin for merging charts into boundary (default: 0.95)
%
% COCO utility functions
%   <a href="matlab: coco_recipes_edit atlas1d_v5 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc atlas1d_v5_demo atlas1d_v5_demo">atlas1d_v5_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
