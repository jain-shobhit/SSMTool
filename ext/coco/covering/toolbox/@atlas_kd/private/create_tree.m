function tree = create_tree(dim)
% CREATE_TREE   Create a tree instance
%
% TREE = CREATE_TREE(DIM)
%
% TREE - Binary tree.
% DIM  - Embedding dimension.
%
% A tree is a structure with a root that is a leaf and a number of charts.
% An octree representation.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

tree.nCharts = 0;                   % Number of charts
tree.dim     = dim;                 % Embedding dimension
tree.root    = create_leaf(dim, 1); % Create instance of leaf

end
