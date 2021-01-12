function A = coco_set_opts(A, B, C)
%COCO_SET_OPTS  Create/alter selected filds in settings structure.
%
%   SETTINGS = COCO_SET_OPTS(SETTINGS, DEFAULTS)
%   SETTINGS = COCO_SET_OPTS(SETTINGS, DEFAULTS, USER)

fprintf(2, '%s: warning: function %s is obsolete, use coco_merge instead.\n', ...
  mfilename, mfilename);

if nargin <= 2
  C = struct();
end

% merge fields from C that are already present in B
fields = fieldnames(B);
for i=1:length(fields)
  field = fields{i};
  
  if isfield(C, field)
    if isempty(B.(field))
      B.(field) = C.(field);
    elseif isstruct(B.(field)) && isstruct(C.(field))
      B.(field) = coco_set(B.(field), C.(field));
    elseif isstruct(B.(field)) && ~isstruct(C.(field))
      error('%s: %s: attempt to overwrite struct with non-struct\n', ...
        mfilename, field);
    elseif ~isstruct(B.(field)) && isstruct(C.(field))
      error('%s: %s: attempt to overwrite non-struct with struct\n', ...
        mfilename, field);
    else
      B.(field) = C.(field);
    end
  end
  
end

% merge new settings into options structure
A = coco_set(A, B);

end
