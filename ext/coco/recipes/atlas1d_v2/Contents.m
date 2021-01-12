% Toolbox: 'atlas1d'
% Version:  2.0
% 
% Source:   Sect. 12.2 of Recipes for Continuation
%
% This toolbox implements a 1-dimensional advancing local cover atlas
% algorithm with varying step size and theta-based projection condition and
% predictor.
%
% Class definition
%   <a href="matlab: coco_recipes_edit atlas1d_v2/@atlas_1d_min atlas_1d_min">atlas_1d_min</a> - 1d atlas class.
%
% Interface methods
%   <a href="matlab: coco_recipes_edit atlas1d_v2/@atlas_1d_min init_prcond">init_prcond</a> - Initialize projection condition.
%   <a href="matlab: coco_recipes_edit atlas1d_v2/@atlas_1d_min init_atlas">init_atlas</a>  - Initialize atlas.
%   <a href="matlab: coco_recipes_edit atlas1d_v2/@atlas_1d_min flush">flush</a>       - Flush point list.
%   <a href="matlab: coco_recipes_edit atlas1d_v2/@atlas_1d_min predict">predict</a>     - Compute predictor.
%   <a href="matlab: coco_recipes_edit atlas1d_v2/@atlas_1d_min add_chart">add_chart</a>   - Add chart to point list.
%   <a href="matlab: coco_recipes_edit atlas1d_v2/@atlas_1d_min refine_step">refine_step</a> - Refine step size.
%
% Private properties
%   base_chart   : Struct used to hold a base chart for the predictor.
%   cont         : Struct used to hold class settings.
%
% Private methods
%   <a href="matlab: coco_recipes_edit atlas1d_v2/@atlas_1d_min get_settings">get_settings</a> - Defines default class settings.
%
% Public methods
%   <a href="matlab:coco_recipes_edit atlas1d_v2/@atlas_1d_min atlas_1d_min>create">create</a> - Static construction method.
%
% Class settings:
% h     : Initial step size (default: 0.1)
% PtMX  : Maximum number of continuation steps (default: 50)
% theta : Theta predictor (default: 0.5)
% hmax  : Maximum step size (default: 0.1)
% hmin  : Minimum step size (default: 0.01)
% hfinc : Factor of step size increment (default: 1.1)
% hfred : Factor of step size decrement (default: 0.5)
% almax : Maximum angle between successive tangent vectors (default: 10)
%
% COCO utility functions
%   <a href="matlab: coco_recipes_edit atlas1d_v2 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc atlas1d_v2_demo atlas1d_v2_demo">atlas1d_v2_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
