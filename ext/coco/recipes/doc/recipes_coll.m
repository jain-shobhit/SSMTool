% The 'coll' toolbox.
% Recipes for Continuation
%
% Constructors and instance functions for trajectory segment continuation
% problem objects associated with autonomous differential equations in
% terms of continuous, piecewise polynomial approximants on a uniform
% discretization of the time domain.
%
% Toolbox development lines:
%   <a href="matlab: coco_recipes_doc coll_v1 coll_v1">coll_v1</a>      - Implementation of collocation method.
%   <a href="matlab: coco_recipes_doc coll_v2 coll_v2">coll_v2</a>      - Supports zero-length initial segments.
%   <a href="matlab: coco_recipes_doc coll_v3 coll_v3">coll_v3</a>      - Error estimation and terminal events.
%   <a href="matlab: coco_recipes_doc coll_v4 coll_v4">coll_v4</a>      - Comoving mesh algorithm with fixed order.
%     +- <a href="matlab: coco_recipes_doc coll_v7 coll_v7">coll_v7</a> - Scaling factor for adaptation step size.
%   <a href="matlab: coco_recipes_doc coll_v5 coll_v5">coll_v5</a>      - Moving mesh algorithm with fixed order.
%   <a href="matlab: coco_recipes_doc coll_v6 coll_v6">coll_v6</a>      - Moving mesh algorithm with varying order.
%
% Toolbox demos:
%   <a href="matlab: coco_recipes_doc coll_v1_demo coll_v1_demo">coll_v1_demo</a> - Command-line use of 'coll' constructors
%   <a href="matlab: coco_recipes_doc coll_v2_demo coll_v2_demo">coll_v2_demo</a> - Initialization with single point solution.
%   <a href="matlab: coco_recipes_doc coll_v3_demo coll_v3_demo">coll_v3_demo</a> - Manual changes to segment discretization.
%   <a href="matlab: coco_recipes_doc coll_v4_demo coll_v4_demo">coll_v4_demo</a> - Error equidistribution with comoving mesh.
%   <a href="matlab: coco_recipes_doc coll_v5_demo coll_v5_demo">coll_v5_demo</a> - Fixed-order adaptive segment discretization.
%   <a href="matlab: coco_recipes_doc coll_v6_demo coll_v6_demo">coll_v6_demo</a> - Variable-order adaptive segment discretization.
%
% See also <a href="matlab: coco_recipes_doc doc recipes_bvp">bvp</a>, <a href="matlab: coco_recipes_doc doc recipes_po">po</a>,  <a href="matlab: coco_recipes_doc doc recipes_msbvp">msbvp</a>, <a href="matlab: coco_recipes_doc doc recipes_hspo">hspo</a>, <a href="matlab: coco_recipes_doc doc recipes_varcoll">varcoll</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: recipes_coll.m 2839 2015-03-05 17:09:01Z fschild $
