function atlas = add_chart2tree(atlas, id)
% ADD_CHART2TREE   Add chart to tree
%
% ATLAS = ADD_CHART2TREE(ATLAS, ID)
%
% ATLAS  - Manifold atlas
% ID     - Chart id (integer)
%
% Adds a chart with a given id to an existing tree.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

[atlas.charts, atlas.tree.root] = ...
  add_chart2leaf(atlas.charts, atlas.tree.root, atlas.tree.dim, id, ...
  atlas.charts{id}.R, atlas.charts{id}.xp);
%atlas.tree.lftBox = atlas.tree.root.lftBox; % checked if used
%atlas.tree.rgtBox = atlas.tree.root.rgtBox; % checked if used

end

function [charts, leaf] = add_chart2leaf(charts, leaf, dim, id, R, center)
% ADD_CHART2LEAF   Add chart to leaf
%
% [ATLAS LEAF] = ADD_CHART2LEAF(ATLAS, LEAF, DIM, ID, R, CENTER)
%
% ATLAS  - Manifold atlas.
% LEAF   - Binary tree leaf.
% DIM    - Manifold dimension (integer).
% ID     - Chart id (integer).
% R      - Radius (real).
% CENTER - Coordinate array (reals).
%
% Adds a chart with a given center and radius to an existing leaf. Function
% is called recursively.

if isempty(leaf.lftTree) && isempty(leaf.rgtTree) && leaf.nCharts<5
  % Insert the new chart, sorted on the d'th coordinate of the center.
  
  % inserts new chart values at beginning of lists
  leaf.values  = [center(leaf.d) leaf.values ];
  leaf.ids     = [      id       leaf.ids    ];
  leaf.radii   = [      R        leaf.radii  ];
  leaf.centers = [   center      leaf.centers];
  
  % sorts all lists based on ascending leaf values
  [leaf.values, idx] = sort(leaf.values, 'ascend');
  leaf.ids           = leaf.ids(idx);
  leaf.radii         = leaf.radii(idx);
  leaf.centers       = leaf.centers(idx);
  
  % increments number of charts
  leaf.nCharts      = leaf.nCharts + 1;
  
  if(leaf.nCharts==1) % No existing bounding box
    leaf.lftBox     = center - R*ones(dim,1);
    leaf.rgtBox     = center + R*ones(dim,1);
  else
    leaf.lftBox     = min(center - R*ones(dim,1), leaf.lftBox);
    leaf.rgtBox     = max(center + R*ones(dim,1), leaf.rgtBox);
  end
  
elseif isempty(leaf.lftTree) && isempty(leaf.rgtTree)
  % Split the leaf and push old and new charts into subtrees
  centers    = cell2mat([center leaf.centers]);
  x1         = max(centers, [], 2);
  x0         = min(centers, [], 2);
  [~, idx]   = sort(x1-x0, 'descend');
  
  leaf.split = (x1(idx(1))+x0(idx(1)))/2;
  leaf.d     = idx(1);
  
  nextdir    = 1;
  if dim>1
    nextdir  = idx(2);
  end
  
  centers     = cell2mat(leaf.centers);
  leaf.values = centers(leaf.d,:);
  
  leaf.lftTree = create_leaf(dim, nextdir);
  leaf.rgtTree = create_leaf(dim, nextdir);
  
  % divide existing leaves into left and right of split
  for i=1:leaf.nCharts
    if leaf.values(i)<leaf.split || (leaf.values(i)==leaf.split && ...
        leaf.lftTree.nCharts<leaf.rgtTree.nCharts)
      charts{leaf.ids(i)}.loc = [charts{leaf.ids(i)}.loc 'l'];
      [charts, leaf.lftTree] = add_chart2leaf(charts, leaf.lftTree, dim, ...
        leaf.ids(i), leaf.radii(i), leaf.centers{i});
    else
      charts{leaf.ids(i)}.loc = [charts{leaf.ids(i)}.loc 'r'];
      [charts, leaf.rgtTree] = add_chart2leaf(charts, leaf.rgtTree, dim, ...
        leaf.ids(i), leaf.radii(i), leaf.centers{i});
    end
  end
  
  if center(leaf.d)<leaf.split || (center(leaf.d)==leaf.split && ...
      leaf.lftTree.nCharts<leaf.rgtTree.nCharts)
    charts{id}.loc = [charts{id}.loc 'l'];
    [charts, leaf.lftTree] = add_chart2leaf(charts, leaf.lftTree, dim, ...
      id, R, center);
  else
    charts{id}.loc = [charts{id}.loc 'r'];
    [charts, leaf.rgtTree] = add_chart2leaf(charts, leaf.rgtTree, dim, ...
      id, R, center);
  end
  leaf.ids     = [];
  leaf.centers = {};
  leaf.values  = [];
  leaf.radii   = [];
  leaf.nCharts = 0;
  
else
  % Push new chart into appropriate subtree.
  if center(leaf.d)<leaf.split || (center(leaf.d)==leaf.split && ...
      leaf.lftTree.nCharts<leaf.rgtTree.nCharts)
    charts{id}.loc = [charts{id}.loc 'l'];
    [charts, leaf.lftTree] = add_chart2leaf(charts, leaf.lftTree, dim, ...
      id, R, center);
  else
    charts{id}.loc = [charts{id}.loc 'r'];
    [charts, leaf.rgtTree] = add_chart2leaf(charts, leaf.rgtTree, dim, ...
      id, R, center);
  end
  
end

% Update bounding box.
if ~isempty(leaf.lftTree)
  leaf.lftBox = leaf.lftTree.lftBox;
  leaf.rgtBox = leaf.lftTree.rgtBox;
end
if ~isempty(leaf.rgtTree)
  leaf.lftBox = min(leaf.lftBox, leaf.rgtTree.lftBox);
  leaf.rgtBox = max(leaf.rgtBox, leaf.rgtTree.rgtBox);
end

end
