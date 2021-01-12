function coco_clear_cache(cmd)
coco_read_solution();
efunc_FDF();
efunc_DFDX();
if nargin>0 && strcmpi(cmd, 'reset')
  coco_func_data.pointers('set', []);
end
end
