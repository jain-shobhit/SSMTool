% Toolbox: 'alg'
% Version:  3.0
% 
% Source:   Sect. 4.1.3 of Recipes for Continuation
%
% This version of the 'alg' toolbox demonstrates the use of a slot function
% to save toolbox data to disk with each labelled solution point.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit alg_v3 alg_construct_eqn">alg_construct_eqn</a> - Construct an instance of 'alg'.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit alg_v3 alg_construct_eqn>alg_F">alg_construct_eqn>alg_F</a>    - COCO-compatible zero function wrapper.
%   <a href="matlab:coco_recipes_edit alg_v3 alg_construct_eqn>alg_DFDU">alg_construct_eqn>alg_DFDU</a> - COCO-compatible linearization of zero
%                                function wrapper.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit alg_v3 alg_read_solution">alg_read_solution</a> - Read 'alg' solution from disk.
%
% Toolbox data (struct)
%   fhan    : Function handle to zero function.
%   dfdxhan : Function handle to Jacobian w.r.t. problem variables.
%   dfdphan : Function handle to Jacobian w.r.t. problem parameters.
%   x_idx   : Index set for problem variables.
%   p_idx   : Index set for problem parameters.
%
% COCO utility functions
%   <a href="matlab: coco_recipes_edit alg_v3 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% see also <a href="matlab: coco_recipes_doc alg_v3_demo alg_v3_demo">alg_v3_demo</a>, coco_save_data

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
