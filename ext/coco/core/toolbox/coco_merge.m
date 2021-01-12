function A = coco_merge(A, B, varargin)
% COCO_MERGE(A, B, [FILTER]) merge (selected) fields of B into A.
% FILTER = cell array of strings
%        | STRING, ...
if isempty(varargin)
  A = coco_opts_tree.merge(A, B);
elseif numel(varargin)==1 && iscell(varargin{1})
  A = coco_opts_tree.merge(A, B, varargin{1});
else
  A = coco_opts_tree.merge(A, B, varargin);
end
end
