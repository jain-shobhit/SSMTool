% Toolbox: 'alg'
% Version:  9.0
% 
% Source:   Sect. 17.1.1 of Recipes for Continuation
%
% This version of the 'alg' toolbox demonstrates the detection of
% bifurcation points using an eigenvalue-based monitor function of type
% regular. Detectable events corresponding to this function and associated
% with 0 are Hopf- and neutral saddle points. A distinction between these
% two types of points is implemented in alg_v10.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit alg_v9 alg_isol2eqn">alg_isol2eqn</a>      - Append 'alg' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_sol2eqn">alg_sol2eqn</a>       - Append 'alg' instance constructed from saved data.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_construct_eqn">alg_construct_eqn</a> - Append an instance of 'alg' to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit alg_v9 alg_construct_eqn>alg_F">alg_construct_eqn>alg_F</a>         - COCO-compatible zero function
%                                     wrapper.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_construct_eqn>alg_DFDU">alg_construct_eqn>alg_DFDU</a>      - COCO-compatible linearization of zero
%                                     function wrapper.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_construct_eqn>alg_fold">alg_construct_eqn>alg_fold</a>      - Fold monitor function.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_construct_eqn>alg_fold_DFDU">alg_construct_eqn>alg_fold_DFDU</a> - Linearization of fold monitor
%                                     function.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_construct_eqn>alg_update">alg_construct_eqn>alg_update</a>    - Update bordering vectors of fold
%                                     monitor function.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_construct_eqn>alg_hopf">alg_construct_eqn>alg_hopf</a>      - Hopf monitor function.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_construct_eqn>alg_bddat">alg_construct_eqn>alg_bddat</a>     - Append to bifurcation data cell array.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit alg_v9 alg_get_settings">alg_get_settings</a>  - Read 'alg' toolbox instance settings.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_arg_check">alg_arg_check</a>     - Basic argument checking for 'alg' toolbox.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_init_data">alg_init_data</a>     - Initialize toolbox data for an instance of 'alg'.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_read_solution">alg_read_solution</a> - Read 'alg' solution and toolbox data from disk.
%   <a href="matlab:coco_recipes_edit alg_v9 alg_fhan_DFDX">alg_fhan_DFDX</a>     - Compute x-derivative of DATA.FHAN(X,P).
%   <a href="matlab:coco_recipes_edit alg_v9 alg_fhan_DFDP">alg_fhan_DFDP</a>     - Compute p-derivative of DATA.FHAN(X,P).
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
%   la_idx1 : Index array used by Hopf monitor function
%   la_idx2 : Index array used by Hopf monitor function
%
% Toolbox settings (in alg field of data)
%   norm    : Boolean flag indicating inclusion of Euclidean norm of
%             problem variables in bifurcation data (default: false).
%   FO      : Flag indicating location and detection of fold points using
%             an embedded ('active') or nonembedded ('regular') monitor
%             function, or no detection (otherwise) (default: 'regular'). 
%   HB      : Boolean flag indicating location and detection of Hopf
%             bifurcation points (default: true).
%
% Continuation parameters
%   test.FO : Track fold monitor function (active or regular).
%   test.HB : Track Hopf monitor function (regular).
%
% Event types
%   FO      : Fold point.
%   HB      : Hopf bifurcation or neutral saddle point.
%
% Solution (struct)
%   x       : Problem variables.
%   p       : Problem parameters.
%   u       : Continuation variables.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit alg_v9 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc alg_v9_demo alg_v9_demo">alg_v9_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
