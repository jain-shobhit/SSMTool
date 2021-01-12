% Toolbox: 'alg'
% Version:  4.0
% 
% Source:   Sect. 4.1.4 of Recipes for Continuation
%
% This version of the 'alg' toolbox demonstrates the implementation of
% generalized toolbox constructors using the <a href="matlab:doc coco_stream">coco_stream</a> utility class for
% flexible argument parsing.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit alg_v4 alg_isol2eqn">alg_isol2eqn</a>      - Construct 'alg' instance from initial data.
%   <a href="matlab:coco_recipes_edit alg_v4 alg_sol2eqn">alg_sol2eqn</a>       - Construct 'alg' instance from saved data.
%   <a href="matlab:coco_recipes_edit alg_v4 alg_construct_eqn">alg_construct_eqn</a> - Construct an instance of 'alg'.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit alg_v4 alg_construct_eqn>alg_F">alg_construct_eqn>alg_F</a>    - COCO-compatible zero function wrapper.
%   <a href="matlab:coco_recipes_edit alg_v4 alg_construct_eqn>alg_DFDU">alg_construct_eqn>alg_DFDU</a> - COCO-compatible linearization of zero
%                                function wrapper.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit alg_v4 alg_arg_check">alg_arg_check</a>     - Basic argument checking for 'alg' toolbox.
%   <a href="matlab:coco_recipes_edit alg_v4 alg_init_data">alg_init_data</a>     - Initialize toolbox data for an instance of 'alg'
%   <a href="matlab:coco_recipes_edit alg_v4 alg_read_solution">alg_read_solution</a> - Read 'alg' solution and toolbox data from disk.
%
% Toolbox data (struct)
%   fhan    : Function handle to zero function.
%   dfdxhan : Empty array or function handle to Jacobian w.r.t. problem
%             variables.
%   dfdphan : Empty array or function handle to Jacobian w.r.t. problem
%             parameters.
%   pnames  : Empty cell array or cell array of string labels for
%             continuation parameters tracking problem parameters.
%   x_idx   : Index set for problem variables.
%   p_idx   : Index set for problem parameters.
%
% Solution (struct)
%   x       : Problem variables.
%   p       : Problem parameters.
%   u       : Continuation variables.
%
% COCO utility functions
%   <a href="matlab: coco_recipes_edit alg_v4 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc alg_v4_demo alg_v4_demo">alg_v4_demo</a>, coco_stream

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
