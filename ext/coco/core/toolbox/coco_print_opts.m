function coco_print_opts(varargin)
%COCO_PRINT_OPTS([fhan] prob [mode] [path])

fhan = 1;
mode = '-inherit';
path = '';
tree = true;

s = coco_stream(varargin{:});
if isnumeric(s.peek)
  fhan = s.get;
end
prob = s.get;
if ischar(s.peek) && coco_opts_tree.is_inherit_mode(s.peek)
  mode = s.get;
end
if ischar(s.peek)
  path = s.get;
  tree = false;
end

if tree
  prob.opts.print_tree(fhan, 'all', true);
else
  assert(any(fhan==[1 2]), ...
    '%s: can display options only to screen', mfilename);
  opts = coco_get(prob, mode, path);
  print_opts_rec(path, opts);
end

end

function print_opts_rec(path, A)

if isempty(path)
  path = 'all';
end
fprintf('%s:\n', path);
disp(A);

if isstruct(A)
  fields = fieldnames(A);
  for i=1:numel(fields)
    field = fields{i};
    if isstruct(A.(field))
      fpath = sprintf('%s.%s', path, field);
      print_opts_rec(fpath, A.(field));
    end
  end
else
end

end
