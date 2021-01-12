function val = coco_bd_val(bd, lab, name)

func = @(x) ~isempty(x) && (x==lab);
lab  = coco_bd_col(bd, 'LAB', false);
row  = cellfun(func, lab);
col  = coco_bd_col(bd, name, false);
val  = col{row};

end
