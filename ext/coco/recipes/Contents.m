% Recipes for Continuation
% Harry Dankowicz and Frank Schilder
%
% * <a href="matlab: coco_recipes_doc doc coco_recipes_lod">Toolboxes and demos ordered by chapter</a>
% * <a href="matlab: coco_recipes_doc doc coco_recipes_dbe">Demos ordered by problem name</a>
% * <a href="matlab: coco_recipes_doc doc coco_recipes_dbf">Demos ordered by figure number</a>
%
% <a href="matlab: coco_recipes_doc doc coco_recipes_howto">How to use this documentation</a>
%
% Tutorial toolboxes and demos:
%   <a href="matlab: coco_recipes_doc doc recipes_alg">alg</a>     - Algebraic equations and basic principles.
%   <a href="matlab: coco_recipes_doc doc recipes_bvp">bvp</a>     - Two-point boundary-value problems.
%   <a href="matlab: coco_recipes_doc doc recipes_coll">coll</a>    - Trajectory collocation with (co)moving meshes.
%   <a href="matlab: coco_recipes_doc doc recipes_compalg">compalg</a> - Coupled systems of algebraic equations.
%   <a href="matlab: coco_recipes_doc doc recipes_dft">dft</a>     - Truncated Fourier series with adaptation.
%   <a href="matlab: coco_recipes_doc doc recipes_hspo">hspo</a>    - Periodic orbits in hybrid dynamical systems.
%   <a href="matlab: coco_recipes_doc doc recipes_msbvp">msbvp</a>   - Multisegment boundary-value problems.
%   <a href="matlab: coco_recipes_doc doc recipes_po">po</a>      - Periodic orbits in smooth dynamical systems.
%   <a href="matlab: coco_recipes_doc doc recipes_varcoll">varcoll</a> - Variational equations and sensitivities.
%
% Tutorial atlas algorithms and demos:
%   <a href="matlab: coco_recipes_doc doc recipes_atlas1d">atlas1d</a> - One-dimensional atlas algorithms.
%   <a href="matlab: coco_recipes_doc doc recipes_atlas2d">atlas2d</a> - Two-dimensional atlas algorithms.
%
% Additional demos:
%   <a href="matlab: coco_recipes_doc cmds_demo cmds_demo">cmds_demo</a>    - COCO's command line interface.
%   <a href="matlab: coco_recipes_doc calcvar_demo calcvar_demo">calcvar_demo</a> - The calculus of variations.
%   <a href="matlab: coco_recipes_doc canard_demo canard_demo">canard_demo</a>  - The Van der Pol canard family.
%   <a href="matlab: coco_recipes_doc cover_demo cover_demo">cover_demo</a>   - Figures for generating cover images.
%
% Recipes utility functions and classes:
%   coco_use_recipes_toolbox  - Add or remove tutorial toolboxes to search path.
%   coco_project_opts_recipes - Defines default settings specific to Recipes for Continuation.
%   coco_recipes_doc          - Show recipes file documentation in help browser.
%   coco_recipes_edit         - Open recipes file in editor.
%   atlas_0d_recipes          - 0d default atlas covering algorithm.
%   atlas_1d_recipes          - 1d default atlas covering algorithm.
%   corr_recipes              - Newton based nonlinear solver.
%   lsol_recipes              - Linear solver.

% COCO core functions:  
%   <a href="matlab: doc coco">coco</a>                -
%   <a href="matlab: doc coco_add_chart_data">coco_add_chart_data</a> -
%   <a href="matlab: doc coco_add_event">coco_add_event</a>      -
%   <a href="matlab: doc coco_add_func">coco_add_func</a>       -
%   <a href="matlab: doc coco_add_func_after">coco_add_func_after</a> -
%   <a href="matlab: doc coco_add_glue">coco_add_glue</a>       -
%   <a href="matlab: doc coco_add_pars">coco_add_pars</a>       -
%   <a href="matlab: doc coco_add_signal">coco_add_signal</a>     -
%   <a href="matlab: doc coco_add_slot">coco_add_slot</a>       -
%   <a href="matlab: doc coco_bd_col">coco_bd_col</a>         -
%   <a href="matlab: doc coco_bd_labs">coco_bd_labs</a>        -
%   <a href="matlab: doc coco_bd_read">coco_bd_read</a>        -
%   <a href="matlab: doc coco_change_func">coco_change_func</a>    -
%   <a href="matlab: doc coco_emit">coco_emit</a>           -
%   <a href="matlab: doc coco_exist">coco_exist</a>          -
%   <a href="matlab: doc coco_ezDFDP">coco_ezDFDP</a>         -
%   <a href="matlab: doc coco_ezDFDX">coco_ezDFDX</a>         -
%   <a href="matlab: doc coco_func_data">coco_func_data</a>      -
%   <a href="matlab: doc coco_get">coco_get</a>            -
%   <a href="matlab: doc coco_get_chart_data">coco_get_chart_data</a> -
%   <a href="matlab: doc coco_get_func_data">coco_get_func_data</a>  -
%   <a href="matlab: doc coco_get_id">coco_get_id</a>         -
%   <a href="matlab: doc coco_merge">coco_merge</a>          -
%   <a href="matlab: doc coco_prob">coco_prob</a>           -
%   <a href="matlab: doc coco_project_opts">coco_project_opts</a>   -
%   <a href="matlab: doc coco_read_solution">coco_read_solution</a>  -
%   <a href="matlab: doc coco_remesh">coco_remesh</a>         -
%   <a href="matlab: doc coco_save_data">coco_save_data</a>      -
%   <a href="matlab: doc coco_set">coco_set</a>            -
%   <a href="matlab: doc coco_set_chart_data">coco_set_chart_data</a> -
%   <a href="matlab: doc coco_set_parival">coco_set_parival</a>    -
%   <a href="matlab: doc coco_stream">coco_stream</a>         -
%   <a href="matlab: doc coco_xchg_pars">coco_xchg_pars</a>      -

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $

% <a href="matlab: coco_recipes_doc doc coco_recipes_iod">Index of toolboxes and examples.</a>
