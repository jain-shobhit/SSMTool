% Toolbox: 'atlas2d'
% Version:  5.0
% 
% Source:   Sect. 14.2 of Recipes for Continuation
%
% This toolbox implements a 2-dimensional expanding
% boundary atlas algorithm with constant step size and theta-based
% projection condition and predictor that recognizes computational domain
% boundaries. The atlas stores a chart network to prevent redundant
% coverage, implements a recursive merger inwards from the boundary, and
% maintains a record of a polygonal boundary, in order to prevent
% premature termination.
%
% Class definition
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min atlas_2d_min">atlas_2d_min</a> - 2d atlas class.
%
% Interface methods
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min init_prcond">init_prcond</a> - Initialize projection condition.
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min init_atlas">init_atlas</a>  - Initialize atlas.
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min flush">flush</a>       - Flush point list.
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min predict">predict</a>     - Compute predictor.
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min add_chart">add_chart</a>   - Add chart to point list.
%
% Private properties
%   boundary     : Integer array to hold ids of boundary charts.
%   charts       : Cell array to hold atlas charts.
%   next_pt      : Chart counter.
%   cont         : Struct used to hold class settings.
%
% Private methods
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min get_settings">get_settings</a> - Defines default class settings.
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min merge">merge</a> - Merge curve segment into atlas.
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min merge_recursive">merge_recursive</a> - Merge chart into chart network.
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min subtract_half_space">subtract_half_space</a> - Subtract half-space from polygon.
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min isneighbor">isneighbor</a> - Check if two charts are neighbors.
%   <a href="matlab: coco_recipes_edit atlas2d_v5/@atlas_2d_min isclose">isclose</a> - Check if two charts are close to each other.
%
% Public methods
%   <a href="matlab:coco_recipes_edit atlas2d_v5/@atlas_2d_min atlas_2d_min>create">create</a> - Static construction method.
%
% Class settings
%   h     : Step size (default: 0.1)
%   PtMX  : Maximum number of continuation steps (default: 50)
%   theta : Theta method (default: 0.5)
%   almax : Critical angle between successive tangent vectors (default: 10)
%   Rmarg : Margin for merging charts into boundary (default: 0.95)
%   Ndirs : Number of available directions (default: 6)
%
% Bifurcation data (identifier: 'atlas')
%   boundary : Integer array to hold ids of boundary charts.
%   charts   : Cell array to hold atlas charts.
%
% Toolbox utility functions
%   <a href="matlab: coco_recipes_edit atlas2d_v5 plot_charts">plot_charts</a>   - Plot atlas as individual charts.
%   <a href="matlab: coco_recipes_edit atlas2d_v5 plot_trisurf">plot_trisurf</a> - Plot atlas as triangulated surface.
%
% COCO utility functions
%   <a href="matlab: coco_recipes_edit atlas2d_v5 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc atlas2d_v5_demo atlas2d_v5_demo">atlas2d_v5_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
