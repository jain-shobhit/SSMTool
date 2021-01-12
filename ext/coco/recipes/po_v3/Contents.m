% Toolbox: 'po'
% Version:  3.0
%
% Compatibility: 'coll', versions 4.0, 5.0, 6.0, and 7.0
% 
% Source:   Sects. 18.3 and 20.2 of Recipes for Continuation
%
% This version of the 'po' toolbox supports adaptive mesh-selection
% algorithms by accommodating the further division of toolbox data
% implemented in the 'coll' toolbox and adaptive changes to the associated
% function dependency index sets during continuation.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit po_v3 po_isol2orb">po_isol2orb</a>  - Append 'po' instance constructed from initial data.
%   <a href="matlab:coco_recipes_edit po_v3 po_sol2orb">po_sol2orb</a>   - Append 'po' instance constructed from saved data.
%   <a href="matlab:coco_recipes_edit po_v3 po_close_orb">po_close_orb</a> - Append an instance of 'po' to problem.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit po_v3 po_close_orb>po_F">po_close_orb>po_F</a>      - Evaluate boundary conditions and
%                            phase condition.
%   <a href="matlab:coco_recipes_edit po_v3 po_close_orb>po_DFDU">po_close_orb>po_DFDU</a>   - Evaluate linearization of
%                            boundary conditions and phase condition.
%   <a href="matlab:coco_recipes_edit po_v3 po_close_orb>po_update">po_close_orb>po_update</a> - Update discretized linear operator.
%   <a href="matlab:coco_recipes_edit po_v3 po_close_orb>po_remesh">po_close_orb>po_remesh</a> - Update dependency index set after remeshing
%                            of 'coll'.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit po_v3 po_init_data">po_init_data</a>     - Initialise toolbox data of po toolbox.
%   <a href="matlab:coco_recipes_edit po_v3 po_read_solution">po_read_solution</a> - Read 'po' solution and toolbox data from disk.
%
% Toolbox data (func_data)
%   x0_idx  : Index array for trajectory end point at t=0.
%   x1_idx  : Index array for trajectory end point at t=1.
%   intfac  : Invariant part of discretized linear operator.
%   xp0     : Discretized linear operator.
%   J       : Linearization of 'po' zero functions.
%   tbid    : Toolbox instance identifier.
%
% Continuation parameters
%   period  : Track interval length (active).
%
% Solution (struct)
%   t       : Array of basepoint time instances.
%   x       : Array of trajectory basepoint values.
%   p       : Problem parameters.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit po_v3 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc coll_v4_demo coll_v4_demo">coll_v4_demo</a>, <a href="matlab: coco_recipes_doc coll_v5_demo coll_v5_demo">coll_v5_demo</a>, <a href="matlab: coco_recipes_doc coll_v6_demo coll_v6_demo">coll_v6_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
