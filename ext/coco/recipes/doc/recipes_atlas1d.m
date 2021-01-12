% The 'atlas1d' toolbox.
% Recipes for Continuation
%
% Class definitions and interface methods for one-dimensional atlas
% algorithms illustrating core principles of atlas expansion and
% consolidation along one-dimensional solution manifolds.
%
% Atlas algorithm development lines:
%   <a href="matlab: coco_recipes_doc atlas1d_v1 atlas1d_v1">atlas1d_v1</a>      - Advancing local cover algorithm.
%     +- <a href="matlab: coco_recipes_doc atlas1d_v2 atlas1d_v2">atlas1d_v2</a> - Step size control and two-step predictor.
%        <a href="matlab: coco_recipes_doc atlas1d_v8 atlas1d_v8">atlas1d_v8</a> - Supports moving mesh adaptation.
%   <a href="matlab: coco_recipes_doc atlas1d_v3 atlas1d_v3">atlas1d_v3</a>      - Expanding boundary algorithm.
%   <a href="matlab: coco_recipes_doc atlas1d_v4 atlas1d_v4">atlas1d_v4</a>      - Terminates at computational boundaries.
%   <a href="matlab: coco_recipes_doc atlas1d_v5 atlas1d_v5">atlas1d_v5</a>      - Can start from computational boundaries.
%   <a href="matlab: coco_recipes_doc atlas1d_v6 atlas1d_v6">atlas1d_v6</a>      - Atlas events using nonembedded monitor functions.
%   <a href="matlab: coco_recipes_doc atlas1d_v7 atlas1d_v7">atlas1d_v7</a>      - Branch switching at branch points.
%
% Atlas algorithm demos:
%   <a href="matlab: coco_recipes_doc atlas1d_v1_demo atlas1d_v1_demo">atlas1d_v1_demo</a> - Regular and singular points and convergence.
%   <a href="matlab: coco_recipes_doc atlas1d_v2_demo atlas1d_v2_demo">atlas1d_v2_demo</a> - Adaptive changes to continuation step size.
%   <a href="matlab: coco_recipes_doc atlas1d_v3_demo atlas1d_v3_demo">atlas1d_v3_demo</a> - Monitoring of atlas boundary and chart gaps.
%   <a href="matlab: coco_recipes_doc atlas1d_v4_demo atlas1d_v4_demo">atlas1d_v4_demo</a> - Interval constraints on the computational domain.
%   <a href="matlab: coco_recipes_doc atlas1d_v5_demo atlas1d_v5_demo">atlas1d_v5_demo</a> - Detection of admissible directions of continuation.
%   <a href="matlab: coco_recipes_doc atlas1d_v6_demo atlas1d_v6_demo">atlas1d_v6_demo</a> - Fold and branch point detection.
%   <a href="matlab: coco_recipes_doc atlas1d_v7_demo atlas1d_v7_demo">atlas1d_v7_demo</a> - Exhaustive branch swithing in nonsmooth oscillator.
%   <a href="matlab: coco_recipes_doc atlas1d_v8_demo atlas1d_v8_demo">atlas1d_v8_demo</a> - Adaptive remeshing of interpolant.
%
% See also <a href="matlab: coco_recipes_doc doc recipes_atlas2d">atlas2d</a>

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: recipes_atlas1d.m 2839 2015-03-05 17:09:01Z fschild $
