function [data res] = coco_save_data(prob, data, varargin) %#ok<INUSL>

if isa(data, 'coco_func_data')
  res = data.data;
else
  res = data;
end

if isfield(res, 'no_save') && iscellstr(res.no_save)
  res = rmfields_rec(res, res.no_save);
end

end

function res = rmfields_rec(res, fields)
for i=1:numel(fields)
  res = rmfield_rec(res, fields{i});
end
end

function res = rmfield_rec(res, field)
[f p] = strtok(field, '.');
if isfield(res, f)
  if isempty(p)
    res = rmfield(res, f);
  else
    res.(f) = rmfield_rec(res.(f), p);
  end
end
end
