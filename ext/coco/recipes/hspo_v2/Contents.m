% Toolbox: 'hspo'
% Version:  2.0
%
% Compatibility: 'msbvp', version 1.0
% 
% Source:   Sects. 17.2.2 and 17.3.3 of Recipes for Continuation
%
% This toolbox implements a multi-segment boundary-value problem for
% periodic solutions of autonomous piecewise smooth ordinary differential
% equations (hybrid dynamical systems). This toolbox demonstrates the
% implementation of monitor functions for bifurcations of hybrid periodic
% solutions as well as the use of event handlers for further processing.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit hspo_v2 hspo_isol2segs">hspo_isol2segs</a> - Append 'hspo' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit hspo_v2 hspo_add_bifus">hspo_add_bifus</a>    - Append monitor functions and events to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit hspo_v2 hspo_bc">hspo_bc</a> - Multi-segment periodic boundary conditions.
%   <a href="matlab:coco_recipes_edit hspo_v2 hspo_bc_DFDX">hspo_bc_DFDX</a> - Linearization of multi-segment periodic boundary conditions.
%   <a href="matlab:coco_recipes_edit hspo_v2 hspo_add_bifus>hspo_TF">hspo_add_bifus>hspo_TF</a>     - Monitor function: bifurcations and stability.
%   <a href="matlab:coco_recipes_edit hspo_v2 hspo_add_bifus>hspo_evhan_NS">hspo_add_bifus>hspo_evhan_NS</a> - Neimark-Sacker bifurcation event handler.
%   <a href="matlab:coco_recipes_edit hspo_v2 hspo_add_bifus>hspo_evhan_PD">hspo_add_bifus>hspo_evhan_PD</a> - Period-doubling bifurcation event handler.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit hspo_v2 hspo_add_bifus>hspo_P">hspo_add_bifus>hspo_P</a> - Compute collection of transfer matrices.
%   <a href="matlab:coco_recipes_edit hspo_v2 hspo_get_settings">hspo_get_settings</a> - Read 'hspo' toolbox instance settings.
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
%   la_idx1 : Index array used by Neimark-Sacker monitor function.
%   la_idx2 : Index array used by Neimark-Sacker monitor function.
%
% Toolbox settings (in hspo field of data)
%   bifus   : Boolean flag indicating detection and location of
%             saddle-node, period-doubling, and Neimark-Sacker bifurcations
%             (default: true).
%   NSad    : Boolean flag indicating location and detection of neutral
%             saddle points (default: false).
%
% Continuation parameters
%   test.SN : Track saddle-node monitor function (regular).
%   test.PD : Track period-doubling monitor function (regular).
%   test.NS : Track Neimark-Sacker monitor function (regular).
%   test.stab : Track number of unstable Floquet multipliers (regular).
%
% Event types
%   SN      : Saddle-node bifurcation.
%   PD      : Period-doubling bifurcation.
%   NS      : Neimark-Sacker bifurcation.
%   NSad    : Neutral saddle point.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit hspo_v2 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc hspo_v2_demo hspo_v2_demo">hspo_v2_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
