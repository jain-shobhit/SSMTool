function varargout = coco_get(prob, varargin)
%COCO_GET(prob, [mode,] [path, [p]... ]) get property P of PATH.

if isempty(prob) || ~isfield(prob, 'opts')
  prob = coco_set(prob);
end

mode = '-inherit';
s = coco_stream(varargin{:});
if coco_opts_tree.is_inherit_mode(s.peek)
  mode = s.get;
end
level = coco_opts_tree.inherit_level(mode);
path  = s.get;
if strcmpi(path, 'all')
  path = '';
end

idx = 1;
while idx==1 || (idx<=nargout && ~isempty(s))
  varargout{idx} = prob.opts.prop_get(path, s.get, level);
  idx = idx+1;
end

end
