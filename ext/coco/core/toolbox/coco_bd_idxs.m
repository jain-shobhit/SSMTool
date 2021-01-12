function idxs = coco_bd_idxs(bd, pt_type)

if nargin<2
  pt_type = [];
end

if ischar(pt_type)
	if strcmpi(pt_type, 'all')
		pt_type = '';
  end
end

if isempty(pt_type)
  labs = coco_bd_col(bd, 'LAB', false);
  f    = @(x) ~isempty(x);
  idxs = find(cellfun(f, labs))';
else
  types = coco_bd_col (bd, 'TYPE');
  idxs  = find(strcmp(pt_type, types))';
end
