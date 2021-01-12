% Toolbox: 'atlas2d'
% Version:  2.0
% 
% Source:   Sect. 13.2 of Recipes for Continuation
%
% This subclass to AtlasBase implements a 2-dimensional point-cloud atlas
% algorithm with constant step size and theta-based projection condition
% and predictor. The atlas stores entire atlas to prevent redundent
% coverage. This algorithm may terminate prematurely.
%
% Class definition
%   <a href="matlab: coco_recipes_edit atlas2d_v2/@atlas_2d_min atlas_2d_min">atlas_2d_min</a> - 2d atlas class.
%
% Interface methods
%   <a href="matlab: coco_recipes_edit atlas2d_v2/@atlas_2d_min init_prcond">init_prcond</a> - Initialize projection condition.
%   <a href="matlab: coco_recipes_edit atlas2d_v2/@atlas_2d_min init_atlas">init_atlas</a>  - Initialize atlas.
%   <a href="matlab: coco_recipes_edit atlas2d_v2/@atlas_2d_min flush">flush</a>       - Flush point list.
%   <a href="matlab: coco_recipes_edit atlas2d_v2/@atlas_2d_min predict">predict</a>     - Compute predictor.
%   <a href="matlab: coco_recipes_edit atlas2d_v2/@atlas_2d_min add_chart">add_chart</a>   - Add chart to point list.
%
% Private properties
%   boundary     : Integer array to hold ids of boundary charts.
%   charts       : Cell array to hold atlas charts.
%   next_pt      : Chart counter.
%   cont         : Struct used to hold class settings.
%
% Private methods
%   <a href="matlab: coco_recipes_edit atlas2d_v2/@atlas_2d_min get_settings">get_settings</a> - Defines default class settings.
%   <a href="matlab: coco_recipes_edit atlas2d_v2/@atlas_2d_min merge">merge</a> - Merge curve segment into atlas.
%   <a href="matlab: coco_recipes_edit atlas2d_v2/@atlas_2d_min isneighbor">isneighbor</a> - Check if two charts are neighbors.
%
% Public methods
%   <a href="matlab:coco_recipes_edit atlas2d_v2/@atlas_2d_min atlas_2d_min>create">create</a> - Static construction method.
%
% Class settings
%   h     : Step size (default: 0.1)
%   PtMX  : Maximum number of continuation steps (default: 50)
%   theta : Theta method (default: 0.5)
%   almax : Critical angle between successive tangent vectors (default: 10)
%   Rmarg : Margin for merging charts into boundary (default: 0.95)
%   Ndirs : Number of available directions (default: 6)
%
% COCO utility functions
%   <a href="matlab: coco_recipes_edit atlas2d_v2 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc atlas2d_v2_demo atlas2d_v2_demo">atlas2d_v2_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
