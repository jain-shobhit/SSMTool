% Toolbox: 'hspo'
% Version:  1.0
%
% Compatibility: 'msbvp', version 1.0
% 
% Source:   Sect. 9.3.1 of Recipes for Continuation
%
% This toolbox demonstrates the implementation of a multi-segment
% boundary-value problem for periodic solutions of autonomous piecewise
% smooth ordinary differential equations (hybrid dynamical systems). The
% toolbox is a wrapper to 'msbvp' and relies on 'msbvp' for a restarter and
% solution extractor. Toolbox data is stored in bc_data field of 'msbvp'
% instance data.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit hspo_v1 hspo_isol2segs">hspo_isol2segs</a> - Append 'hspo' instance constructed from initial data.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit hspo_v1 hspo_bc">hspo_bc</a> - Multi-segment periodic boundary conditions.
%   <a href="matlab:coco_recipes_edit hspo_v1 hspo_bc_DFDX">hspo_bc_DFDX</a> - Linearization of multi-segment periodic boundary conditions.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit hspo_v1 hspo_arg_check">hspo_arg_check</a> - Basic argument checking for hspo toolbox.
%   <a href="matlab:coco_recipes_edit hspo_v1 hspo_init_data">hspo_init_data</a> - Initialise toolbox data of hspo toolbox.
%
% Toolbox data (struct or func_data [contingent on bc_update])
%   hhan    : Function handle to event function.
%   dhdxhan : Optional function handle to Jacobian of event function w.r.t.
%             problem variables.
%   dhdphan : Optional function handle to Jacobian of event function w.r.t.
%             problem parameters.
%   ghan    : Function handle to reset function.
%   dgdxhan : Optional function handle to Jacobian of reset function w.r.t.
%             problem variables.
%   dgdphan : Optional function handle to Jacobian of reset function w.r.t.
%             problem parameters.
%   modes   : Sequence of mode identifiers.
%   events  : Sequence of event identifiers.
%   resets  : Sequence of reset identifiers.
%   nsegs   : Number of segments.
%   x0_idx  : Index arrays for trajectory end points at t=0.
%   x1_idx  : Index arrays for trajectory end points at t=1.
%   dim     : Array of state-space dimensions.
%   cdim    : Total state-space dimensions.
%   pdim    : Number of problem parameters.
%   rows    : Row index array for sparse Jacobian.
%   cols    : Column index array for sparse Jacobian.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit hspo_v1 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc hspo_v1_demo hspo_v1_demo">hspo_v1_demo</a>, <a href="matlab:coco_recipes_doc msbvp_v1 msbvp_v1">msbvp_v1</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
