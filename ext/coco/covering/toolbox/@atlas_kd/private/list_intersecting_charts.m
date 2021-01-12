function list = list_intersecting_charts(tree, id, center, R)
%LIST_INTERSECTING_CHARTS   Finds the charts which overlap the given center and radius.
%
% LIST = LIST_INTERSECTING_CHARTS(TREE, CHART, CENTER, R)
%
% LIST   - Array of chart id's (integers).
% TREE   - Binary tree.
% CHART  - Chart id (integer).
% CENTER - Point (Array of reals).
% R      - Radius (real).
%
% In the tree given by the tree input argument, find all charts, except
% that with id given by the chart input argument, that overlap a circle
% with radius R centered at the point whose coordinates are given by the
% center input argument.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

list = get_intersecting_charts(tree.root, tree.dim, [], id, center, R);
end

function list = get_intersecting_charts(leaf, dim, list, id, center, R)
%GET_INTERSECTING_CHARTS   Append charts that overlap in a leaf.
%
% LIST = GET_INTERSECTING_CHARTS(LEAF, K, LIST, CHART, CENTER, R)
%
% Append to the list input argument all charts in the leaf given by the
% leaf input argument, except that with id given by the chart input
% argument, that overlap a circle with radius R centered at the point whose
% coordinates are given by the center input argument. Function is called
% recursively.

if isempty(leaf) % safe exit
  return
  
elseif isempty(leaf.lftTree) && leaf.nCharts==0 % reached an empty leaf
  return
  
elseif any(center - R*ones(dim,1) - leaf.rgtBox>0) || ... % no overlap
    any(center + R*ones(dim,1) - leaf.lftBox<0)
  return
  
elseif isempty(leaf.lftTree) % no branches
  for i=1:leaf.nCharts
    if leaf.ids(i)~=id
      if norm(center - leaf.centers{i}) < R + leaf.radii(i)
        list = [list leaf.ids(i)]; %#ok<AGROW>
      end
    end
  end
  
else % otherwise check left and right subtrees
  list = get_intersecting_charts(leaf.lftTree, dim, list, id, center, R);
  list = get_intersecting_charts(leaf.rgtTree, dim, list, id, center, R);
end

end
