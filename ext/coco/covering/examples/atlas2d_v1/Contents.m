% Toolbox: 'atlas2d'
% Version:  1.0
% 
% Source:   Sect. 13.1.1 of Recipes for Continuation
%
% This toolbox implements a 2-dimensional point-cloud atlas algorithm with
% constant step size and theta-based projection condition and predictor.
% This basic algorithm typically double covers and may terminate
% prematurely. The algorithm generalizes the source in Recipes for
% Continuation by providing support for a projected geometry onto a subset
% of elements of the vector of continuation variables.
%
% Class definition
%   <a href="matlab: edit @atlas_2d_min/atlas_2d_min">atlas_2d_min</a> - 2d atlas class.
%
% Interface methods
%   <a href="matlab: edit @atlas_2d_min/init_prcond">init_prcond</a> - Initialize projection condition.
%   <a href="matlab: edit @atlas_2d_min/init_atlas">init_atlas</a>  - Initialize atlas.
%   <a href="matlab: edit @atlas_2d_min/flush">flush</a>       - Flush point list.
%   <a href="matlab: edit @atlas_2d_min/predict">predict</a>     - Compute predictor.
%   <a href="matlab: edit @atlas_2d_min/add_chart">add_chart</a>   - Add chart to point list.
%
% Private properties
%   boundary     : Cell array to hold atlas boundary charts.
%   next_pt      : Chart counter.
%   cont         : Struct used to hold class settings.
%
% Private methods
%   <a href="matlab: edit @atlas_2d_min/get_settings">get_settings</a> - Defines default class settings.
%   <a href="matlab: edit @atlas_2d_min/merge">merge</a> - Merge curve segment into atlas.
%   <a href="matlab: edit @atlas_2d_min/isneighbor">isneighbor</a> - Check if two charts are neighbors.
%
% Public methods
%   <a href="matlab: edit @atlas_2d_min/atlas_2d_min>create">create</a> - Static construction method.
%
% Class settings
%   h     : Step size (default: 0.1)
%   PtMX  : Maximum number of continuation steps (default: 50)
%   theta : Theta method (default: 0.5)
%   almax : Critical angle between successive tangent vectors (default: 10)
%   Rmarg : Margin for merging charts into boundary (default: 0.95)
%   Ndirs : Number of available directions (default: 6)
%
% See also <a href="matlab: doc ../atlas2d_v1_demo/Contents.m">atlas2d_v1_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 3116 2019-06-09 02:11:55Z hdankowicz $
