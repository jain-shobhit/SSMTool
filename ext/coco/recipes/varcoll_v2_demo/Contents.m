% Toolbox: 'varcoll'
% Version:  2.0
% 
% Source:   Sect. 10.2.2 of Recipes for Continuation
%
% The example in varcoll_v2_demo illustrates the use of an embedded
% variational zero problem and multiple instances of 'coll', with
% appropriate boundary conditions, to continue approximate heteroclinic
% connecting orbits between an equilibrium and a periodic orbit.
%
% Scripts
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo demo_lorenz">demo_lorenz.m</a>       - Figures 10.3, 10.4, and 10.5.
%
% Functions
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo figure_10_1">figure_10_1.m</a>       - Generate figure 10.1.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo figure_10_2">figure_10_2.m</a>       - Generate figure 10.2.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo figure_10_3">figure_10_3.m</a>       - Generate figure 10.3.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo figure_10_4">figure_10_4.m</a>       - Generate figure 10.4.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo figure_10_5">figure_10_5.m</a>       - Generate figure 10.5.
%
% Other files
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo lorenz">lorenz.m</a>            - 'coll'-compatible encoding of vector field.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo lorenz_DFDX">lorenz_DFDX.m</a>       - 'coll'-compatible encoding of Jacobian w.r.t. states.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo lorenz_DFDP">lorenz_DFDP.m</a>       - 'coll'-compatible encoding of Jacobian w.r.t. problem parameters.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo lorenz_DFDXDX">lorenz_DFDXDX.m</a>     - 'varcoll'-compatible encoding of second partials w.r.t. states.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo lorenz_DFDXDP">lorenz_DFDXDP.m</a>     - 'varcoll'-compatible encoding of mixed second partials.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo povar_isol2orb">povar_isol2orb.m</a>    - Append 'po' and 'varcoll' instances from initial data.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo povar_sol2orb">povar_sol2orb.m</a>     - Append 'po' and 'varcoll' instance from stored data.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo riess_start_1">riess_start_1.m</a>     - Append additional constrained 'coll' instances from initial data.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo riess_restart_1">riess_restart_1.m</a>   - Append additional constrained 'coll' instances from stored data.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo riess_close_het_1">riess_close_het_1.m</a> - Append eigenspace and boundary conditions to 'po' and 'coll' instances.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo riess_start_2">riess_start_2.m</a>     - Append heteroclinic connection problem from initial data.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo riess_restart_2">riess_restart_2.m</a>   - Append heteroclinic connection problem from stored data.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo riess_close_het_2">riess_close_het_2.m</a> - Append additional boundary conditions to constrained 'coll' instances.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo eig_bcs">eig_bcs.m</a>           - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo lingap">lingap.m</a>            - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo linphase">linphase.m</a>          - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo proj_bcs">proj_bcs.m</a>          - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit varcoll_v2_demo var_evs">var_evs.m</a>           - COCO-compatible function encoding.
%
% See also <a href="matlab: coco_recipes_doc varcoll_v2 varcoll_v2">varcoll_v2</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
