% Toolbox: 'atlas2d'
% Version:  5.0
% 
% Source:   Sect. 14.2 of Recipes for Continuation
%
% This toolbox implements a 2-dimensional expanding boundary atlas
% algorithm with constant step size and theta-based projection condition
% and predictor that recognizes computational domain boundaries. The atlas
% stores a chart network to prevent redundant coverage, implements a
% recursive merger inwards from the boundary, and maintains a record of a
% polygonal boundary, in order to prevent premature termination. The
% algorithm generalizes the source in Recipes for Continuation by providing
% support for a projected geometry onto a subset of elements of the vector
% of continuation variables.
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
%   boundary     : Integer array to hold ids of boundary charts.
%   charts       : Cell array to hold atlas charts.
%   next_pt      : Chart counter.
%   cont         : Struct used to hold class settings.
%
% Private methods
%   <a href="matlab: edit @atlas_2d_min/get_settings">get_settings</a> - Defines default class settings.
%   <a href="matlab: edit @atlas_2d_min/merge">merge</a> - Merge curve segment into atlas.
%   <a href="matlab: edit @atlas_2d_min/merge_recursive">merge_recursive</a> - Merge chart into chart network.
%   <a href="matlab: edit @atlas_2d_min/subtract_half_space">subtract_half_space</a> - Subtract half-space from polygon.
%   <a href="matlab: edit @atlas_2d_min/isneighbor">isneighbor</a> - Check if two charts are neighbors.
%   <a href="matlab: edit @atlas_2d_min/isclose">isclose</a> - Check if two charts are close to each other.
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
% Bifurcation data (identifier: 'atlas')
%   boundary : Integer array to hold ids of boundary charts.
%   charts   : Cell array to hold atlas charts.
%
% Toolbox utility functions
%   <a href="matlab: edit plot_charts">plot_charts.m</a>   - Plot atlas as individual charts.
%   <a href="matlab: edit plot_trisurf">plot_trisurf.m</a> - Plot atlas as triangulated surface.
%
% See also <a href="matlab: doc ../atlas2d_v5_demo/Contents.m">atlas2d_v5_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 3123 2019-06-11 03:14:29Z hdankowicz $
