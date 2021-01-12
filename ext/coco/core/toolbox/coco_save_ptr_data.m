function [data_ptr res] = coco_save_ptr_data(opts, data_ptr, varargin) %#ok<INUSL>

fprintf(2, '%s: slot function will become obsolete, use ''coco_save_data'' instead\n', ...
  mfilename);

res = data_ptr.data;

if isfield(res, 'no_save')
  res = rmfield(res, res.no_save);
end
