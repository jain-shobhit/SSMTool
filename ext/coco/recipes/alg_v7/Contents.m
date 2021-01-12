% Toolbox: 'alg'
% Version:  7.0
% 
% Source:   Sect. 16.2.2 of Recipes for Continuation
%
% This version of the 'alg' toolbox demonstrates the fundamentals of
% detection and continuation of codimension-1 bifurcation points with the
% example of fold detection and continuation.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit alg_v7 alg_isol2eqn">alg_isol2eqn</a>      - Append 'alg' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_sol2eqn">alg_sol2eqn</a>       - Append 'alg' instance constructed from saved data.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_construct_eqn">alg_construct_eqn</a> - Append an instance of 'alg' to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit alg_v7 alg_construct_eqn>alg_F">alg_construct_eqn>alg_F</a>         - COCO-compatible zero function
%                                     wrapper.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_construct_eqn>alg_DFDU">alg_construct_eqn>alg_DFDU</a>      - COCO-compatible linearization of zero
%                                     function wrapper.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_construct_eqn>alg_fold">alg_construct_eqn>alg_fold</a>      - Fold monitor function.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_construct_eqn>alg_fold_DFDU">alg_construct_eqn>alg_fold_DFDU</a> - Linearization of fold monitor
%                                     function.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_construct_eqn>alg_update">alg_construct_eqn>alg_update</a>    - Update bordering vectors of fold
%                                     monitor function.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_construct_eqn>alg_bddat">alg_construct_eqn>alg_bddat</a>     - Append to bifurcation data cell array.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit alg_v7 alg_get_settings">alg_get_settings</a>  - Read 'alg' toolbox instance settings.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_arg_check">alg_arg_check</a>     - Basic argument checking for 'alg' toolbox.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_init_data">alg_init_data</a>     - Initialize toolbox data for an instance of 'alg'.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_read_solution">alg_read_solution</a> - Read 'alg' solution and toolbox data from disk.
%   <a href="matlab:coco_recipes_edit alg_v7 alg_fhan_DFDX">alg_fhan_DFDX</a>     - Compute x-derivative of DATA.FHAN(X,P).
%   <a href="matlab:coco_recipes_edit alg_v7 alg_fhan_DFDP">alg_fhan_DFDP</a>     - Compute p-derivative of DATA.FHAN(X,P).
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
%   b       : Bordering vector close to left singular vector of Jacobian
%             w.r.t. problem variables with smallest singular value.
%   c       : Bordering vector close to right singular vector of Jacobian
%             w.r.t. problem variables with smallest singular value.
%   rhs     : Right-hand side matrix for bordering algorithm.
%
% Toolbox settings (in alg field of data)
%   norm    : Boolean flag indicating inclusion of Euclidean norm of
%             problem variables in bifurcation data (default: false).
%   FO      : Flag indicating location and detection of fold points using
%             an embedded ('active') or nonembedded ('regular') monitor
%             function, or no detection (otherwise) (default: 'regular'). 
%
% Continuation parameters
%   test.FO : Track fold monitor function.
%
% Event types
%   FO      : Fold point.
%
% Solution (struct)
%   x       : Problem variables.
%   p       : Problem parameters.
%   u       : Continuation variables.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit alg_v7 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc alg_v7_demo alg_v7_demo">alg_v7_demo</a>, coco_xchg_pars

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
