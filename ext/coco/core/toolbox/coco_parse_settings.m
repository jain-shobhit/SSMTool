function sets = coco_parse_settings(prob, spec, curr, tbid)
flds = spec(:,1);
type = spec(:,2);
defs = spec(:,3);
acts = spec(:,4);
args = spec(:,5);

defaults = struct();
for i=1:numel(flds)
  defaults.(flds{i}) = defs{i};
end

sets = coco_get(prob, tbid);
sets = coco_merge(curr, sets);
sets = coco_merge(defaults, sets);

for i=1:numel(flds)
  fld = flds{i};
  typ = type{i};
  act = acts{i};
  arg = args{i};
  
  val = sets.(fld);
  
  switch typ
    case 'log|on'
      flag = (islogical(val) || any(strcmpi(val, {'on' 'off'})));
      assert(flag, '%s: input to setting ''%s'' must be logical or ''on''/''off''', ...
        mfilename, fld);
      
    case '@'
      flag = isa(val, 'function_handle');
      assert(flag, '%s: input to setting ''%s'' must be a function handle', ...
        mfilename, fld);
      
    case '@|[]'
      flag = (isempty(val) || isa(val, 'function_handle'));
      assert(flag, '%s: input to setting ''%s'' must be a function handle or empty', ...
        mfilename, fld);
      
    case 'int'
      flag = (isnumeric(val) && isscalar(val) && mod(val,1)==0);
      assert(flag, '%s: input to setting ''%s'' must be an integer number', ...
        mfilename, fld);
      
    case '[int]'
      flag = (isnumeric(val) && all(mod(val,1)==0));
      assert(flag, '%s: input to setting ''%s'' must be a vector of integer numbers', ...
        mfilename, fld);
      
    case 'num'
      flag = (isnumeric(val) && isscalar(val));
      assert(flag, '%s: input to setting ''%s'' must be a number', ...
        mfilename, fld);
      
    case '[num]'
      flag = isnumeric(val);
      assert(flag, '%s: input to setting ''%s'' must be a numerical array', ...
        mfilename, fld);
      
    case 'str'
      flag = ischar(val);
      assert(flag, '%s: input to setting ''%s'' must be a string', ...
        mfilename, fld);
      
    case '{str}'
      flag = iscellstr(val);
      assert(flag, '%s: input to setting ''%s'' must be a cell array of strings', ...
        mfilename, fld);

    otherwise
      error('%s: illegal type specification ''%s''', mfilename, typ);
  end
  
  if ischar(act)
    
    switch lower(act)
      
      case 'read'
        
      case 'switch'
        if ischar(val)
          flags = strcmpi(val, {'on' 'off'});
          sets.(fld) = flags(1);
        end
        
      case 'numel'
        chk_numel(fld, typ, val, arg{:});
        
      otherwise
        error('%s: invalid action ''%''', mfilename, act);
    end
    
  elseif isa(act, 'function_handle')
    sets.(fld) = act(fld, sets.(fld), defs{i}, arg{:});
    
  else
    error('%s: invalid action specification', mfilename);
  end
end

end

function chk_numel(fld, typ, val, varargin)
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
    msg = sprintf('vector of %s integer numbers', msg);
  case '[num]'
    msg = sprintf('numerical array with %s elements', msg);
  case '{str}'
    msg = sprintf('cell array of %s strings', msg);
  otherwise
    msg1 = sprintf('cannot apply test ''numel'' to setting ''%s''', fld);
    msg2 = sprintf('numel is applicable to types ''[int]'', ''[num]'' and ''{str}''');
    msg3 = sprintf('type of setting %s is ''%s''', fld, typ);
    error('%s: %s\n%s\n%s', mfilename, msg1, msg2, msg3);
end
assert(any(numel(val)==nums), ...
  '%s: input to setting ''%s'' must be a %s', mfilename, fld, msg);
end
