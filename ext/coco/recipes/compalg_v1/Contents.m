% Toolbox: 'compalg'
% Version:  1.0
% 
% Source:   Sects. 5.2.1-5.2.2 of Recipes for Continuation
%
% This version of the 'compalg' toolbox demonstrates embedding of several
% 'alg' instances into a single instance of a composite continuation
% problem.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit compalg_v1 compalg_isol2sys">compalg_isol2sys</a>      - Append and glue several 'alg' instances 
%                           constructed from initial data.
%   <a href="matlab:coco_recipes_edit compalg_v1 compalg_sol2sys">compalg_sol2sys</a>       - Append and glue several 'alg' instances
%                           constructed from saved data.
%   <a href="matlab:coco_recipes_edit compalg_v1 compalg_close_sys">compalg_close_sys</a>     - Add gluing conditions to embedded instances
%                           of 'alg'.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit compalg_v1 compalg_arg_check">compalg_arg_check</a>     - Basic argument checking for 'compalg' toolbox.
%   <a href="matlab:coco_recipes_edit compalg_v1 compalg_read_solution">compalg_read_solution</a> - Read 'compalg' solution and toolbox data
%                           from disk.
%
% Toolbox data (struct)
%   neqs     : Number of 'alg' instances.
%   fhans    : Cell array of function handles to zero functions.
%   dfdxhans : Empty cell array or cell array of empty arrays or function
%              handles to Jacobians w.r.t. problem variables.
%   dfdphan  : Empty cell array or cell array of empty arrays or function
%              handles to Jacobians w.r.t. problem parameters.
%   pnames   : Empty cell array or cell array of string labels for
%              continuation parameters tracking problem parameters.
%
% Solution (struct)
%   x        : Cell array of problem variables.
%   p        : Problem parameters.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit compalg_v1 coco_rm_this_path">coco_rm_this_path</a>     - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc compalg_v1_demo compalg_v1_demo">compalg_v1_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
