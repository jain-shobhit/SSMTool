function [leaf, found] = change_radius(leaf, dim, chart, level, newR, found)
% CHANGE_RADIUS Change the radius of a given chart
%
% [LEAF FOUND] = CHANGE_RADIUS(LEAF, DIM, CHART, LEVEL, NEWR, FOUND)
%
%   LEAF  - Current subtree
%   DIM   - Embedding dimension
%   CHART - Chart object
%   LEVEL - Current tree depth
%   newR  - New radius
%   FOUND - Boolean variable that indicates location
%
% Function is called recursively to locate chart in tree structure (usually
% called initially with root). Once found, radius is changed and the
% bounding box is modified at all dependent levels.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

at_level = false;

if level<=numel(chart.loc) && ~found
  if chart.loc(level) == 'l'
    [leaf.lftTree, found] = change_radius(leaf.lftTree, dim, chart, ...
      level+1, newR, found);
  elseif chart.loc(level) == 'r'
    [leaf.rgtTree, found] = change_radius(leaf.rgtTree, dim, chart, ...
      level+1, newR, found);
  end
elseif level>numel(chart.loc) && ~found
  idx = find(leaf.ids==chart.id);
  if ~isempty(idx)
    leaf.radii(idx) = newR;
    at_level        = true;
    found           = true;
  end
end

if found && at_level
  lftBox = cell2mat(leaf.centers) - repmat(leaf.radii, dim, 1);
  leaf.lftBox = min(lftBox, [], 2);
  rgtBox = cell2mat(leaf.centers) + repmat(leaf.radii, dim, 1);
  leaf.rgtBox = max(rgtBox, [], 2);
elseif found && ~at_level
  leaf.lftBox = min(leaf.lftTree.lftBox, leaf.rgtTree.lftBox);
  leaf.rgtBox = max(leaf.lftTree.rgtBox, leaf.rgtTree.rgtBox);
end

end
