% Toolbox: 'alg'
% Version:  6.0
% 
% Source:   Sect. 5.3 of Recipes for Continuation
%
% This version of the 'alg' toolbox demonstrates the support of toolbox
% settings and implements a user-defined slot function.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit alg_v6 alg_isol2eqn">alg_isol2eqn</a>      - Append 'alg' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit alg_v6 alg_sol2eqn">alg_sol2eqn</a>       - Append 'alg' instance constructed from saved data.
%   <a href="matlab:coco_recipes_edit alg_v6 alg_construct_eqn">alg_construct_eqn</a> - Append an instance of 'alg' to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit alg_v6 alg_construct_eqn>alg_F">alg_construct_eqn>alg_F</a>     - COCO-compatible zero function wrapper.
%   <a href="matlab:coco_recipes_edit alg_v6 alg_construct_eqn>alg_DFDU">alg_construct_eqn>alg_DFDU</a>  - COCO-compatible linearization of zero
%                                function wrapper.
%   <a href="matlab:coco_recipes_edit alg_v6 alg_construct_eqn>alg_bddat">alg_construct_eqn>alg_bddat</a> - Append to bifurcation data cell array.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit alg_v6 alg_get_settings">alg_get_settings</a>  - Read 'alg' toolbox instance settings.
%   <a href="matlab:coco_recipes_edit alg_v6 alg_arg_check">alg_arg_check</a>     - Basic argument checking for 'alg' toolbox.
%   <a href="matlab:coco_recipes_edit alg_v6 alg_init_data">alg_init_data</a>     - Initialize toolbox data for an instance of 'alg'.
%   <a href="matlab:coco_recipes_edit alg_v6 alg_read_solution">alg_read_solution</a> - Read 'alg' solution and toolbox data from disk.
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
% Toolbox settings (in alg field of data)
%   norm    : Boolean flag indicating inclusion of Euclidean norm of
%             problem variables in bifurcation data (default: false).
%
% Solution (struct)
%   x       : Problem variables.
%   p       : Problem parameters.
%   u       : Continuation variables.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit alg_v6 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc alg_v6_demo alg_v6_demo">alg_v6_demo</a>, coco_set, coco_get, coco_merge

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
