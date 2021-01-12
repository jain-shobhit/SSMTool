function val = coco_bd_vals(bd, labs, name)

func = @(x) ~isempty(x) && any(x==labs);
labs = coco_bd_col(bd, 'LAB', false);
row  = cellfun(func, labs);
col  = coco_bd_col(bd, name, false);
val  = col(row);

try
  val2 = val;
  if size(val2{1},1)>1
    for i=1:numel(val2)
      val2{i} = reshape(val2{i}, [1 size(val2{i})]);
    end
  end
  val = cat(1, val2{:});
catch e %#ok<NASGU>
end

end
