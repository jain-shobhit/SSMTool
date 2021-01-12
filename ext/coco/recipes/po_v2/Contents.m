% Toolbox: 'po'
% Version:  2.0
%
% Compatibility: 'coll', versions 1.0 and 2.0
% 
% Source:   Sect. 17.2.1 and 17.3.2 of Recipes for Continuation
%
% This version of the 'po' toolbox demonstrates the implementation of
% monitor functions for bifurcations of periodic solutions as well as the
% use of event handlers for further processing.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit po_v2 po_isol2orb">po_isol2orb</a>  - Append 'po' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit po_v2 po_sol2orb">po_sol2orb</a>   - Append 'po' instance constructed from saved data.
%   <a href="matlab:coco_recipes_edit po_v2 po_close_orb">po_close_orb</a> - Append one instance of po to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit po_v2 po_close_orb>po_F')">po_close_orb>po_F</a>      - Evaluate boundary conditions and
%                            phase condition.
%   <a href="matlab:coco_recipes_edit po_v2 po_close_orb>po_DFDU">po_close_orb>po_DFDU</a>   - Evaluate linearization of
%                            boundary conditions and phase condition.
%   <a href="matlab:coco_recipes_edit po_v2 po_close_orb>po_update">po_close_orb>po_update</a> - Update discretized linear operator.
%   <a href="matlab:coco_recipes_edit po_v2 po_close_orb>po_TF">po_close_orb>po_TF</a>     - Update discretized linear operator.
%   <a href="matlab:coco_recipes_edit po_v2 po_close_orb>po_evhan_NS">po_close_orb>po_evhan_NS</a> - Neimark-Sacker bifurcation event handler.
%   <a href="matlab:coco_recipes_edit po_v2 po_close_orb>po_evhan_PD">po_close_orb>po_evhan_PD</a> - Period-doubling bifurcation event handler.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit po_v2 po_get_settings">po_get_settings</a>  - Read 'po' toolbox instance settings.
%   <a href="matlab:coco_recipes_edit po_v2 po_init_data">po_init_data</a>     - Initialize toolbox data for an instance of 'po'.
%   <a href="matlab:coco_recipes_edit po_v2 po_read_solution">po_read_solution</a> - Read 'po' solution and toolbox data from disk.
%
% Toolbox data (func_data)
%   x0_idx  : Index array for trajectory end point at t=0.
%   x1_idx  : Index array for trajectory end point at t=1.
%   intfac  : Invariant part of discretized linear operator.
%   xp0     : Discretized linear operator.
%   J       : Linearization of 'po' zero functions.
%   la_idx1 : Index array used by Neimark-Sacker monitor function.
%   la_idx2 : Index array used by Neimark-Sacker monitor function.
%   tbid    : Toolbox instance identifier.
%
% Toolbox settings (in po field of data)
%   bifus   : Boolean flag indicating detection and location of
%             saddle-node, period-doubling, and Neimark-Sacker bifurcations
%             (default: true).
%   NSad    : Boolean flag indicating location and detection of neutral
%             saddle points (default: false).
%
% Continuation parameters
%   period  : Track interval length (active).
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
% Solution (struct)
%   t       : Array of basepoint time instances.
%   x       : Array of trajectory basepoint values.
%   p       : Problem parameters.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit po_v2 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc po_v2_demo po_v2_demo">po_v2_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
