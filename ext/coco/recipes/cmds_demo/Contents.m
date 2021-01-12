% Toolbox: 'coco' 
%
% Source:   Sects. 3.3 and 15.3 of Recipes for Continuation
%
% The examples in cmds_demo illustrate the basic command-line interface for
% constructing a continuation problem structure and performing continuation
% of the solutions to a suitable restricted continuation problem. The
% examples demonstrate the use of the <a href="matlab: doc coco_add_func">coco_add_func</a> constructor, the
% special-purpose wrappers <a href="matlab: doc coco_add_pars">coco_add_pars</a> and <a href="matlab: doc coco_add_glue">coco_add_glue</a>, as well as the
% application of the <a href="matlab: doc coco">coco</a> entry-point function, and the <a href="matlab: doc coco_get_func_data">coco_get_func_data</a>,
% <a href="matlab: doc coco_add_event">coco_add_event</a>, <a href="matlab: doc coco_read_solution">coco_read_solution</a>, and <a href="matlab: doc coco_bd_read">coco_bd_read</a> utilities.
%
% Scripts
%   <a href="matlab: coco_recipes_edit cmds_demo demo_basic">demo_basic.m</a>       - Examples 3.7, 3.9, and 3.10.
%   <a href="matlab: coco_recipes_edit cmds_demo demo_stages">demo_stages.m</a>      - Examples 3.11 and 3.12.
%   <a href="matlab: coco_recipes_edit cmds_demo demo_henon">demo_henon.m</a>       - Example 3.13.
%   <a href="matlab: coco_recipes_edit cmds_demo demo_brusselator">demo_brusselator.m</a> - Examples 3.14 and 4.9.
%   <a href="matlab: coco_recipes_edit cmds_demo demo_events">demo_events.m</a>      - Example 15.1.
%
% Functions
%   <a href="matlab: coco_recipes_edit cmds_demo figure_3_2">figure_3_2.m</a>       - Generate figure 3.2.
%   <a href="matlab: coco_recipes_edit cmds_demo figure_3_3">figure_3_3.m</a>       - Generate figure 3.3.
%
% Other files
%   <a href="matlab: coco_recipes_edit cmds_demo circ">circ.m</a>             - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit cmds_demo dist">dist.m</a>             - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit cmds_demo plan">plan.m</a>             - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit cmds_demo hype">hype.m</a>             - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit cmds_demo henon">henon.m</a>            - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit cmds_demo period_A">period_A.m</a>         - Period-n orbits, minimal version.
%   <a href="matlab: coco_recipes_edit cmds_demo period_B">period_B.m</a>         - Period-n orbits, redundant variables.
%   <a href="matlab: coco_recipes_edit cmds_demo finitediff">finitediff.m</a>       - COCO-compatible function encoding.
%   <a href="matlab: coco_recipes_edit cmds_demo euclid">euclid.m</a>           - COCO-compatible function encoding.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: Contents.m 2839 2015-03-05 17:09:01Z fschild $
