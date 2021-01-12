function leaf = create_leaf(dim, d)
% CREATE_LEAF   Creates a leaf in a tree
%
% LEAF = CREATE_LEAF(DIM, D)
%
% LEAF - Binary tree leaf.
% DIM  - Embedding dimension.
% D    - Coordinate direction of split.
%
% A leaf is a subtree. The topmost tree is the root of a tree structure.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

leaf.d       = d;                  % Coordinate direction of split
leaf.split   = 0.;                 % Coordinate value of split
leaf.nCharts = 0;                  % Number of charts at top level of leaf

leaf.ids     = [];                 % Array of chart ids at top level of leaf
leaf.radii   = [];                 % Array of chart radii
leaf.centers = cell(0);            % Cell array of chart centers

leaf.values  = [];                 % Array of split values of each chart

leaf.lftBox  = 1.e20*ones(dim,1);  % Lower left corner of this leaf's box
leaf.rgtBox  = -1.e20*ones(dim,1); % Upper right corner of the box
leaf.lftTree = [];                 % Left subtree
leaf.rgtTree = [];                 % Right subtree

end
