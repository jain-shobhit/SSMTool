function opts = coco_set_func_data(opts, fid, varargin)

error('%s: obsolete function', mfilename);

efunc = opts.efunc;
if any(strcmp(fid, efunc.identifyers))
  for funcarray = { 'zero' 'embedded' 'regular' 'singular' }
    for i=1:numel(efunc.(funcarray{1}))
      func = efunc.(funcarray{1})(i);
      if strcmp(fid, func.identifyer)
        iarg = 1;
        while iarg<=nargin-2
          switch lower(varargin{iarg})
            case 'data'
              func.data  = varargin{iarg+1};
            case 'cdata'
              func.cdata = varargin{iarg+1};
            case 'xidx'
              func.x_idx = varargin{iarg+1};
            case 'fidx'
              func.f_idx = varargin{iarg+1};
            case 'midx'
              func.m_idx = varargin{iarg+1};
            case 'x0'
              opts.efunc.x0(func.x_idx) = varargin{iarg+1};
            otherwise
              error('%s: unknown function data field ''%s''', ...
                mfilename, varargin{iarg});
          end
          iarg = iarg + 2;
        end
        efunc.(funcarray{1})(i) = func;
        opts.efunc              = efunc;
        return
      end
    end
  end
else
  error('%s: could not find function with identifyer ''%s''', mfilename, fid);
end
