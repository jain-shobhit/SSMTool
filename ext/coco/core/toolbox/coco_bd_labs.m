function labs = coco_bd_labs(bd, pt_type)

if nargin<2
  pt_type = [];
end

if ischar(pt_type)
	if strcmpi(pt_type, 'all')
		pt_type = '';
  end
end

labs = coco_bd_col(bd, 'LAB', false);

if isempty(pt_type)
  labs = [ labs{:} ];
else
  types = coco_bd_col (bd, 'TYPE');
  idx   = strcmp(pt_type, types);
  labs  = [ labs{idx} ];
end
