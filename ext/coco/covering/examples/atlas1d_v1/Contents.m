% Toolbox: 'atlas1d'
% Version:  1.0
% 
% Source:   Sect. 12.1.1 of Recipes for Continuation
%
% This toolbox implements a 1-dimensional advancing local cover atlas
% algorithm with constant step size and pseudo-arclength-based projection
% condition and predictor. The algorithm generalizes the source in Recipes
% for Continuation by providing support for a projected geometry onto a
% subset of elements of the vector of continuation variables.
%
% Class definition
%   <a href="matlab: edit @atlas_1d_min/atlas_1d_min">atlas_1d_min</a> - 1d atlas class.
%
% Interface methods
%   <a href="matlab: edit @atlas_1d_min/init_prcond">init_prcond</a>  - Initialize projection condition.
%   <a href="matlab: edit @atlas_1d_min/init_atlas">init_atlas</a>   - Initialize atlas.
%   <a href="matlab: edit @atlas_1d_min/flush">flush</a>        - Flush point list.
%   <a href="matlab: edit @atlas_1d_min/predict">predict</a>      - Compute predictor.
%   <a href="matlab: edit @atlas_1d_min/add_chart">add_chart</a>    - Add chart to point list.
%
% Private properties
%   base_chart   : Struct used to hold a base chart for the predictor.
%   cont         : Struct used to hold class settings.
%
% Private methods
%   <a href="matlab: edit @atlas_1d_min/get_settings">get_settings</a> - Defines default class settings.
%
% Public methods
%   <a href="matlab: edit @atlas_1d_min/atlas_1d_min>create">create</a> - Static construction method.
%
% Class settings
%   h            : Step size (default: 0.1)
%   PtMX         : Maximum number of continuation steps (default: 50)
%
% See also <a href="matlab: doc ../atlas1d_v1_demo/Contents.m">atlas1d_v1_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 3110 2019-06-08 04:25:59Z hdankowicz $
