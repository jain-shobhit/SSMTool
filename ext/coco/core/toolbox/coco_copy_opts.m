function subsettings = coco_copy_opts(settings, varargin)

fprintf(2, '%s: warning: function %s is obsolete, use coco_merge instead.\n', ...
  mfilename, mfilename);

if isstruct(varargin{1})
  % varargin = { new_settings fieldname ... }
  %          | { new_settings fieldlist }
  % copy of selected fields
  new_settings = varargin{1};
  if ischar(varargin{2})
    fields = varargin(2:end);
  else
    fields = varargin{2};
  end
  subsettings = settings;
  for i=1:numel(fields)
    field = fields{i};
    if isfield(new_settings, field)
      subsettings.(field) = new_settings.(field);
    end
  end
else
  % varargin = { fieldname ... } | { fieldlist }
  % create structure with copy of selected fields
  if ischar(varargin{1})
    fields = varargin;
  else
    fields = varargin{1};
  end
  subsettings = struct();
  for i=1:numel(fields)
    field = fields{i};
    if isfield(settings, field)
      subsettings.(field) = settings.(field);
    end
  end
end
