function data = coco_parse_data(grammar, spec, varargin)
%COCO_PARSE_DATA   Construct data from input according to grammar.
%
% This is a very primitive parser. To enhance functionality, add a compile
% function that compiles the grammar into a proper parsing tree
% specification. The curent algorithm is sufficient for most uses.

toks = spec(:,1);
peek = spec(:,2);
type = spec(:,3);
flds = spec(:,4);
defs = spec(:,5);
acts = spec(:,6);
args = spec(:,7);

expand = {};
data = struct();
for i=1:numel(flds)
  fld = flds{i};
  if ~isempty(fld)
    dval = defs{i};
    n = ischar(dval)*numel(dval);
    if n>=1 && dval(1)=='$'
      if n>=2 && dval(2)=='$'
        data.(fld) = dval(2:end);
      else
        data.(fld) = dval;
        expand = [ expand fld ]; %#ok<AGROW>
      end
    else
      data.(fld) = dval;
    end
  end
end

str  = coco_stream(varargin{:});
token(grammar);
data = parse_rec(data, toks, peek, type, flds, defs, acts, args, str);

for i=1:numel(expand)
  fld = expand{i};
  val = data.(fld);
  n = ischar(val)*numel(val);
  if n>=1 && val(1)=='$'
    j = find(strcmp(val(2:end), toks), 1, 'first');
    data.(fld) = data.(flds{j});
  end
end
end

function data = parse_rec(data, toks, peek, type, flds, defs, acts, args, str)
while true
  [tok idx flag] = read_tok_check_val(toks, peek, type, str);
  switch lower(tok)
    case '['
      [tok idx flag] = read_tok_check_val(toks, peek, type, str);
      if isempty(idx)
        switch lower(tok)
          case 'opts'
            skip_rec();
            return;
          otherwise
            error('%s: unsupported grammar specification ''... [ %s ...', ...
              mfilename, tok);
        end
      elseif flag
        data.(flds{idx}) = str.get(peek{idx});
        data = parse_rec(data, toks, peek, type, flds, defs, acts, args, str);
      else
        skip_rec();
      end
      
    case ']'
      return;
      
    case { '{' '}' '|' '(' ')' }
      error('%s: parse operator ''%s'' not implemented', mfilename, tok);
      
    case { '' 'opts' }
      return;
      
    otherwise
      if isempty(idx)
        error('%s: unsupported grammar specification', mfilename);
      elseif ~flag
        val_error(idx, tok, type, peek);
      end
      act = acts{idx};
      fld = flds{idx};
    
      if ischar(act)
        
        switch lower(act)
          
          case 'read'
            data.(fld) = str.get(peek{idx});
            
          otherwise
            error('%s: invalid action ''%''', mfilename, act);
        end
        
      elseif isa(act, 'function_handle')
        
        arg = args{idx};
        data.(fld) = act(tok, opts.(fld), defs{idx}, arg{:});
        
      else
        error('%s: invalid action specification', mfilename);
      end
  end
end

end

function skip_rec()
while true
  tok = token();
  switch tok
    case ']'
      return
    case '['
      skip_rec();
    case ''
      error('%s: invalid grammar string', mfilename);
  end
end
end

function [tok idx flag] = read_tok_check_val(toks, peek, type, str)
tok = token();
idx = find(strcmpi(tok, toks), 1, 'first');
if isempty(idx) || strcmpi(tok, 'opts')
  flag = false;
  return;
end
val = str.peek(peek{idx});
switch(type{idx})
  
  case '@'
    flag = isa(val, 'function_handle');
    
  case '{@}'
    flag = (iscell(val) && ...
      all(cellfun(@(x) isa(x, 'function_handle'), val)));
    
  case '@|[]'
    flag = (isempty(val) || isa(val, 'function_handle'));
    
  case '{@|[]}'
    flag = (iscell(val) && ...
      all(cellfun(@(x) (isempty(x) || isa(x, 'function_handle')), val)));
    
  case 'num'
    flag = (isnumeric(val) && isscalar(val));
    
  case '[num]'
    flag = isnumeric(val);
    
  case 'num|@f'
    flag = ((isnumeric(val) && isscalar(val)) || ...
      isa(val, 'function_handle'));
    
  case '{[num]}'
    flag = iscell(val) && all(cellfun(@(x) isnumeric(x), val));
    
  case 'str'
    flag = ischar(val);
    
  case '{str}'
    flag = iscellstr(val);
    
  case 'stct'
    flag = isstruct(val) || isa(val,'coco_func_data');
    
  otherwise
    error('%s: illegal type specification ''%s'' in grammar string', mfilename, type{idx});
end
end

function val_error(idx, tok, type, peek)
switch(type{idx})
  
  case '@'
    msg = 'not a function handle';
    
  case '{@}'
    msg = 'not a cell array of function handles';
    
  case '@|[]'
    msg = 'not a function handle or empty';
    
  case '{@|[]}'
    msg = 'not a cell array of function handles or empties';
    
  case 'num'
    msg = 'not a number';
    
  case '[num]'
    msg = 'not a numerical array';
    
  case 'num|@f'
    msg = 'not a number or function handle';
    
  case '{[num]}'
    msg = 'not a cell array of numerical arrays';
    
  case 'str'
    msg = 'not a string';
    
  case '{str}'
    if strcmpi(peek{idx}, 'cell')
      msg = 'not a string or cell array of strings';
    else
      msg = 'not a cell array of strings';
    end
    
  case 'stct'
    msg = 'not a struct or coco_func_data object';
    
  otherwise
    error('%s: illegal type specification ''%s'' in grammar string', mfilename, type{idx});
end
error('%s: argument for token ''%s'' %s', mfilename, tok, msg);
end

function tok = token(new_grammar)
persistent grammar pos
if nargin>0
  grammar = [ new_grammar ' '];
  pos     = 1;
  return;
end
while pos<=numel(grammar) && grammar(pos)==' '
  pos = pos+1;
end
if pos>numel(grammar)
  tok = '';
  return;
end
switch grammar(pos)
  case { '[' ']' '|' '{' '}' '(' ')' }
    tok = grammar(pos);
    pos = pos+1;
  otherwise
    tok = '';
    while ~any(grammar(pos)==' []|{}()')
      tok = [ tok grammar(pos) ]; %#ok<AGROW>
      pos = pos+1;
    end
end
end
