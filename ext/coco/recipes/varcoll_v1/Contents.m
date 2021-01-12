% Toolbox: 'varcoll'
% Version:  1.0
% 
% Source:   Sect. 10.1 of Recipes for Continuation
%
% This toolbox implements a monitor function that computes the fundamental
% solution to the first variational equation along a segment object. This
% toolbox demonstrates attaching functions to existing instances of other
% independently developed toolboxes.
%
% The use of the var_coll_add toolbox constructor in Section 10.1 of
% Recipes for Continuation is poorly motivated and the associated
% implementation of the po_mult_eigs_bddat and hspo_mult_eigs_bddat slot
% functions is incorrect, since it makes use of function data that may not
% be associated with the chart being flushed to disk. The two alternative
% implementations in this folder correct this error.
%
% Toolbox constructors
%   <a href="matlab:coco_recipes_edit varcoll_v1 var_coll_add">var_coll_add</a>  - Attach 'varcoll' instance to segment object.
%   <a href="matlab:coco_recipes_edit varcoll_v1 po_mult_add">po_mult_add</a>   - Add slot and/or monitor function computing Floquet
%                   multipliers of single-segment periodic orbit.
%   <a href="matlab:coco_recipes_edit varcoll_v1 hspo_mult_add">hspo_mult_add</a> - Add slot and/or monitor function computing Floquet
%                   multipliers of hybrid periodic orbit.
%
% Toolbox instance functions
%   <a href="matlab:coco_recipes_edit varcoll_v1 var_coll_add>var_coll_seg">var_coll_add>var_coll_seg</a>          - Compute fundamental solution of
%                                        variational problem.
%   <a href="matlab:coco_recipes_edit varcoll_v1 po_mult_add>po_mult_eigs_bddat">po_mult_add>po_mult_eigs_bddat</a>     - Add Floquet multipliers to
%                                        bifurcation data. Implementation
%                                        in Recipes for Continuation is
%                                        incorrect.
%   <a href="matlab:coco_recipes_edit varcoll_v1 po_mult_add>po_mult_eigs">po_mult_add>po_mult_eigs</a>           - Add Floquet multipliers to
%                                        nonembedded continuation parameters.
%   <a href="matlab:coco_recipes_edit varcoll_v1 hspo_mult_add>hspo_mult_eigs_bddat">hspo_mult_add>hspo_mult_eigs_bddat</a> - Add Floquet multipliers to
%                                        bifurcation data. Implementation
%                                        in Recipes for Continuation is
%                                        incorrect.
%   <a href="matlab:coco_recipes_edit varcoll_v1 hspo_mult_add>hspo_mult_eigs">hspo_mult_add>hspo_mult_eigs</a>       - Add Floquet multipliers to
%                                        nonembedded continuation parameters.
%
% Toolbox utility functions
%   <a href="matlab:coco_recipes_edit varcoll_v1 var_coll_add>var_coll_init_data">var_coll_add>var_coll_init_data</a> - Initialize toolbox data for an
%                                     instance of 'varcoll'.
%   <a href="matlab:coco_recipes_edit varcoll_v1 hspo_mult_add>hspo_P">hspo_mult_add>hspo_P</a>            - Compute collection of transfer
%                                     matrices.
%   <a href="matlab:coco_recipes_edit varcoll_v1 var_eval_sol">var_eval_sol</a>                    - Compute Jacobian of time-T flow.
%
% Toolbox data (struct)
%   tbid   : 'coll' object instance identifier.
%   dim    : State-space dimension.
%   M1_idx : Index array for Jacobian of time-T flow.
%   row    : Initial condition and right-hand side.
%   M      : Fundamental solution to variational problem.
%
% COCO utility functions
%   <a href="matlab:coco_recipes_edit varcoll_v1 coco_rm_this_path">coco_rm_this_path</a> - Remove toolbox from search path.
%
% See also <a href="matlab: coco_recipes_doc varcoll_v1_demo varcoll_v1_demo">varcoll_v1_demo</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
