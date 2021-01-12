function opts = coco_parse_opts(spec, varargin)
opts = struct();
if isempty(spec)
  return;
end

toks = spec(:,1);
flds = spec(:,2);
defs = spec(:,3);
acts = spec(:,4);
args = spec(:,5);

for i=1:numel(flds)
  fld = flds{i};
  if ~isempty(fld)
    opts.(fld) = defs{i};
  end
end

str = coco_stream(varargin{:});
while ~isempty(str)
  tok = str.peek;
  if ischar(tok)
    i = find(strcmpi(tok, toks), 1, 'first');
    if isempty(i)
      break;
    end
    act = acts{i};
    fld = flds{i};
    
    if ischar(act)
      
      switch lower(act)
        
        case { 'break' 'end' }
          str.skip;
          break;
          
        case 'toggle'
          str.skip;
          opts.(fld) = ~defs{i};
          
        case 'switch'
          str.skip;
          val = str.get;
          if ischar(val)
            flags = strcmpi(val, {'on' 'off'});
            assert(any(flags), ...
              '%s: argument to switch ''%s'' must be ''on'' or ''off''', ...
              mfilename, tok);
            val = flags(1);
          elseif islogical(val)
          else
            error('%s: argument to switch ''%s'' must be logical or ''on''/''off''', ...
              mfilename, tok);
          end
          opts.(fld) = val;
          
        case 'read'
          str.skip;
          opts.(fld) = str.get;
          
        otherwise
          error('%s: invalid action ''%''', mfilename, act);
      end
      
    elseif isa(act, 'function_handle')
      arg = args{i};
      opts.(fld) = act(tok, opts.(fld), defs{i}, arg{:});
      
    else
      error('%s: invalid action specification', mfilename);
    end
    
  else
    break;
  end
end

end
