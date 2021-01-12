% The 'atlas2d' toolbox.
% Recipes for Continuation
%
% Class definitions and interface methods for two-dimensional atlas
% algorithms illustrating core principles of atlas expansion and
% consolidation along two-dimensional solution manifolds.
%
% Atlas algorithm development lines:
%   <a href="matlab: coco_recipes_doc atlas2d_v1 atlas2d_v1">atlas2d_v1</a> - Point cloud algorithm.
%   <a href="matlab: coco_recipes_doc atlas2d_v2 atlas2d_v2">atlas2d_v2</a> - Chart network.
%   <a href="matlab: coco_recipes_doc atlas2d_v3 atlas2d_v3">atlas2d_v3</a> - Recursive merge algorithm.
%   <a href="matlab: coco_recipes_doc atlas2d_v4 atlas2d_v4">atlas2d_v4</a> - Henderson's algorithm.
%   <a href="matlab: coco_recipes_doc atlas2d_v5 atlas2d_v5">atlas2d_v5</a> - Terminates at computational boundaries.
%   <a href="matlab: coco_recipes_doc atlas2d_v6 atlas2d_v6">atlas2d_v6</a> - Can start from computational boundaries.
%
% Atlas algorithm demos:
%   <a href="matlab: coco_recipes_doc atlas2d_v1_demo atlas2d_v1_demo">atlas2d_v1_demo</a> - Redundant coverage and premature termination.
%   <a href="matlab: coco_recipes_doc atlas2d_v2_demo atlas2d_v2_demo">atlas2d_v2_demo</a> - Premature termination but no redundant coverage.
%   <a href="matlab: coco_recipes_doc atlas2d_v3_demo atlas2d_v3_demo">atlas2d_v3_demo</a> - Recursive network traversal during consolidation.
%   <a href="matlab: coco_recipes_doc atlas2d_v4_demo atlas2d_v4_demo">atlas2d_v4_demo</a> - Construction of polygonal atlas boundary.
%   <a href="matlab: coco_recipes_doc atlas2d_v5_demo atlas2d_v5_demo">atlas2d_v5_demo</a> - Interval constraints on the computational domain.
%   <a href="matlab: coco_recipes_doc atlas2d_v6_demo atlas2d_v6_demo">atlas2d_v6_demo</a> - Detection of admissible directions of continuation.
%
% See also <a href="matlab: coco_recipes_doc doc recipes_atlas1d">atlas1d</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: recipes_atlas2d.m 2839 2015-03-05 17:09:01Z fschild $
