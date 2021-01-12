classdef coco_stream < handle
  % COCO_STREAM convert cell array to stream object.
  
  properties ( Access = private )
    tokens = {}
    ntok   = 0;
    pos    = 1;
  end
  
  methods
    
    function p = coco_stream(varargin)
      if nargin>=1 && isa(varargin{1}, 'coco_stream')
        p        = varargin{1};
        p.tokens = [ p.tokens varargin{2:end} ];
      else
        p.tokens = varargin;
      end
      p.ntok = numel(p.tokens);
    end
    
    function t = peek(p, varargin)
      if p.pos<=p.ntok
        t = p.tokens{p.pos};
      else
        t = [];
      end
      if nargin>=2 && ~iscell(t) && strcmp('cell', varargin{1})
        t = { t };
      end
    end
    
    function varargout = get(p, varargin)
      for i=1:max(1,nargout)
        varargout{i} = p.peek(varargin{:});
        p.pos        = min(p.pos,p.ntok)+1;
      end
    end
    
    function p = put(p, varargin)
      p.tokens = [p.tokens(1:(p.pos-1)) varargin p.tokens(p.pos:end)];
      p.ntok   = numel(p.tokens);
    end
    
    function skip(p, n)
      if nargin>1
        p.pos = min(p.pos+n,p.ntok+1);
      else
        p.pos = min(p.pos+1,p.ntok+1);
      end
    end
    
    function flag = isempty(p)
      flag = p.pos>p.ntok;
    end
    
    function n = numel(p)
      n = p.ntok-p.pos+1;
    end
    
    function idx = tell(p)
      idx = p.pos;
    end
    
    function p = seek(p, idx)
      p.pos = max(1, min(idx, p.ntok+1));
    end
    
    function disp(p)
      disp('stream =');
      disp(p.tokens(p.pos:end))
    end
  end
  
end
