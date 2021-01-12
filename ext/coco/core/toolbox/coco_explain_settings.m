%function sets = coco_explain_settings(prob, spec, curr, tbid)
function coco_explain_settings(spec, sets)
flds = spec(:,1);
type = spec(:,2);
defs = spec(:,3);
acts = spec(:,4);
args = spec(:,5);
desc = spec(:,6);

width = 5;

defaults = struct();
for i=1:numel(flds)
  fld = flds{i};
  defaults.(fld) = defs{i};
  width = min(15, max(width, length(fld)));
end

% sets = coco_get(prob, tbid);
% sets = coco_merge(curr, sets);
% sets = coco_merge(defaults, sets);

for i=1:numel(flds)
  fld = flds{i};
  typ = type{i};
  act = acts{i};
  arg = args{i};
  des = desc{i};
  
  val = sets.(fld);
  v = 'INVALID INPUT';
  
  switch typ
    
    case 'log|on'
      flag = (islogical(val) || any(strcmpi(val, {'on' 'off'})));
      if flag
        if islogical(val)
          if val
            v = 'true';
          else
            v = 'false';
          end
        else
          v = [ '''' val '''' ];
        end
      end
      t = 'logical or ''on''/''off''';
      fprintf('%*s : %s, %s, {%s}\n', width, fld, des, t, v);
      
    case '@'
      flag = isa(val, 'function_handle');
      if flag
        v = ['@' func2str(val)];
      end
      t = 'function handle';
      fprintf('%*s : %s, %s, {%s}\n', width, fld, des, t, v);
      
    case '@|[]'
      flag = (isempty(val) || isa(val, 'function_handle'));
      if flag
        if isempty(val)
          v = '[]';
        else
          v = ['@' func2str(val)];
        end
      end
      t = 'function handle or empty';
      fprintf('%*s : %s, %s, {%s}\n', width, fld, des, t, v);
      
    case 'int'
      flag = (isnumeric(val) && isscalar(val) && mod(val,1)==0);
      if flag
        v = sprintf('%d', val);
      end
      t = 'integer number';
      fprintf('%*s : %s, %s, {%s}\n', width, fld, des, t, v);
      
    case '[int]'
      flag = (isnumeric(val) && all(mod(val,1)==0));
      t = 'vector of integer numbers';
      if flag
        fprintf('%*s : %s, %s\n', width, fld, des, t);
      else
        fprintf('%*s : %s, %s, {%s}\n', width, fld, des, t, v);
      end
      
    case 'num'
      flag = (isnumeric(val) && isscalar(val));
      if flag
        v = sprintf('%.6g', val);
      end
      t = 'numeric (a number)';
      fprintf('%*s : %s, %s, {%s}\n', width, fld, des, t, v);
      
    case '[num]'
      flag = isnumeric(val);
      t = 'numerical array';
      if flag
        fprintf('%*s : %s, %s\n', width, fld, des, t);
      else
        fprintf('%*s : %s, %s, {%s}\n', width, fld, des, t, v);
      end
      
    case 'str'
      flag = ischar(val);
      if flag
        v = val;
      end
      t = 'character string';
      fprintf('%*s : %s, %s, {%s}\n', width, fld, des, t, v);
      
    case '{str}'
      flag = iscellstr(val);
      t = 'cell array of strings';
      if flag
        fprintf('%*s : %s, %s\n', width, fld, des, t);
      else
        fprintf('%*s : %s, %s, {%s}\n', width, fld, des, t, v);
      end
      
    otherwise
      error('%s: illegal type specification ''%s''', mfilename, typ);
  end
  
  if ischar(act)
    
    switch lower(act)
      
      case 'read'
        
      case 'switch'
        
      case 'numel'
        chk_numel(width, fld, typ, val, arg{:});
        
      otherwise
        error('%s: invalid action ''%''', mfilename, act);
    end
    
  elseif isa(act, 'function_handle')
    fprintf('%*s   customized action @%s\n', width, '', func2str(act));
    
  else
    error('%s: invalid action specification', mfilename);
  end
end

end

function chk_numel(width, fld, typ, val, varargin)
nums = [ varargin{:} ];
assert(~isempty(nums), ...
  '%s: argument for action ''numel'' must not be empty', ...
  mfilename);
msg = sprintf('%d', nums(1));
if numel(nums)>1
  for i=2:numel(nums)-1
    msg = sprintf('%s, %d', msg, nums(i));
  end
  msg = sprintf('%s or %d', msg, nums(end));
end
switch typ
  case '[int]'
    fprintf('%*s   accepting vector of %s integer numbers, ', width, '', msg);
    fprintf('{[');
    if numel(val)>=1
      fprintf(' %d', val(1));
      for i=2:min(numel(val), 4)
        fprintf(', %d', val(i));
      end
      if numel(val)>i
        fprintf(', ... ');
      else
        fprintf(' ');
      end
    end
    fprintf(']}\n');
  case '[num]'
    fprintf('%*s   accepting numerical array with %s elements, ', width, '', msg);
    fprintf('{[');
    if numel(val)>=1
      fprintf(' %.4g', val(1));
      for i=2:min(numel(val), 4)
        fprintf(', %.4g', val(i));
      end
      if numel(val)>i
        fprintf(', ... ');
      else
        fprintf(' ');
      end
    end
    fprintf(']}\n');
  case '{str}'
    fprintf('%*s   accepting cell array of %s strings\n', width, '', msg);
    fprintf('%*s   {{', width, '');
    if numel(val)>=1
      len = 8+length(val{1});
      fprintf(' ''%s''', val{1});
      prdots = false;
      for i=2:min(numel(val), 4)
        len = len+4+length(val{i});
        if len>70
          prdots = true;
          break;
        end
        fprintf(', ''%s''', val{i});
      end
      if prdots || numel(val)>i
        fprintf(', ... ');
      else
        fprintf(' ');
      end
    end
    fprintf('}}\n');
  otherwise
    msg1 = sprintf('cannot apply test ''numel'' to setting ''%s''', fld);
    msg2 = sprintf('numel is applicable to types ''[int]'', ''[num]'' and ''{str}''');
    msg3 = sprintf('type of setting %s is ''%s''', fld, typ);
    error('%s: %s\n%s\n%s', mfilename, msg1, msg2, msg3);
end
assert(any(numel(val)==nums), ...
  '%s: input to setting ''%s'' must be a %s', mfilename, fld, msg);
end
