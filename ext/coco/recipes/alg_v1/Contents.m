% Toolbox: 'alg'
% Version:  1.0
% 
% Source:   Sect. 4.1.1 of Recipes for Continuation
%
% This version of the 'alg' toolbox demonstrates the implementation of a
% basic problem constructor and the use of a zero function wrapper for
% user-friendly definition of a zero problem.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit alg_v1 alg_construct_eqn">alg_construct_eqn</a> - Construct an instance of 'alg'.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit alg_v1 alg_construct_eqn>alg_F">alg_construct_eqn>alg_F</a> - COCO-compatible zero function wrapper.
%
% Toolbox data (struct)
%   x_idx : Index set for problem variables.
%   p_idx : Index set for problem parameters.
%
% COCO utility functions
%   <a href="matlab: coco_recipes_edit alg_v1 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc alg_v1_demo alg_v1_demo">alg_v1_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
