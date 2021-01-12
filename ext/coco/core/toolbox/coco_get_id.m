function id = coco_get_id(prefix, varargin)
%COCO_GET_FID   Create identifyer from prefix and namespace(s).

if isempty(prefix)
  id = '';
  format = '%s%s';
else
  if ~ischar(prefix)
    error('%s: argument prefix must be a string', mfilename);
  end
  id = prefix;
  format = '%s.%s';
end

if numel(varargin)==1 && iscell(varargin{1})
  list = varargin{1};
  ids = cell(size(list));
  for i=1:numel(list)
    name = list{i};
    if isempty(name)
      ids{i} = id;
    else
      if ~ischar(name)
        error('%s: namespace arguments must be strings', mfilename);
      end
      ids{i} = sprintf(format, id, name);
    end
  end
  id = ids;
else
  for i=1:numel(varargin)
    name = varargin{i};
    if ~isempty(name)
      if ~ischar(name)
        error('%s: namespace arguments must be strings', mfilename);
      end
      id = sprintf(format, id, name);
      format = '%s.%s';
    end
  end
end

end
