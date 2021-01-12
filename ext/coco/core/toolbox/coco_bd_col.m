function [col] = coco_bd_col(bd, name, cat_flag)

if nargin<3
  cat_flag = true;
end

if iscell(name)
  for i=1:numel(name)
    col{i} = coco_get_col(bd, name{i}, cat_flag); %#ok<AGROW>
  end
  if cat_flag
    col = cat_arrays(col);
  end
else
  col = coco_get_col(bd, name, cat_flag);
end

end

function col = coco_get_col(bd, name, cat_flag)
col = strcmp(name, { bd{1,:} });
if ~any(col)
  error('%s: column not found: ''%s''', mfilename, name);
end

col = bd(2:end, find(col,1));
if isempty(col)
  return
end

if cat_flag
  col = merge_arrays(col);
end
end

function col = merge_arrays(col)
if all(cellfun(@isnumeric, col)) || all(cellfun(@isempty, col))
  cat_flag = true;
  for d = 1:ndims(col{1})
    if any(cellfun('size', col, d) ~= size(col{1},d))
      cat_flag = false;
      break
    end
  end
  if cat_flag
    if d==2
      [m n] = size(col{1});
      if n == 1
        col = horzcat(col{:});
      elseif m == 1
        col = vertcat(col{:});
      else
        col = cat(d+1, col{:});
      end
    else
      col = cat(d+1, col{:});
    end
  end
end
end

function col = cat_arrays(col)
if all(cellfun(@isnumeric, col)) || all(cellfun(@isempty, col))
  cat_flag = true;
  for d = 2:ndims(col{1})
    if any(cellfun('size', col, d) ~= size(col{1},d))
      cat_flag = false;
      break
    end
  end
  if cat_flag
    col = cat(1, col{:});
  end
end
end
