function varargout = coco_get_func_data(prob, fid, varargin)

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_get_func_data.m 3177 2020-02-20 17:04:54Z hdankowicz $

efunc     = prob.efunc;
varargout = {};
idx       = find(strcmpi(fid, ['efunc' efunc.identifyers]), 1)-1;
if isempty(idx)
  error('%s: could not find function with identifyer ''%s''', mfilename, fid);
else
  if idx
    func = efunc.funcs(idx);
    for oarg = 1:nargin-2
      switch lower(varargin{oarg})
        case 'cdata' % chart data - not in use according to Harry
          varargout = [ varargout { func.data.cdata } ]; %#ok<AGROW>
        case 'data'  % function data structure
          if isa(func.data, 'coco_func_data')
            varargout = [ varargout { func.data.protect() } ]; %#ok<AGROW>
          else
            varargout = [ varargout { func.data } ]; %#ok<AGROW>
          end
        case {'xidx' 'uidx'} % position of x0 and t0 in full solution vector
          varargout = [ varargout { func.x_idx(:) } ]; %#ok<AGROW>
        case 'type' % function type
          varargout = [ varargout { func.type } ]; %#ok<AGROW>
        case 'fidx'  % position of F(x) in zero function
          varargout = [ varargout { func.f_idx(:) } ]; %#ok<AGROW>
        case 'midx'  % position of F(x) in monitor function
          varargout = [ varargout { func.m_idx(:) } ]; %#ok<AGROW>
        case {'x0' 'u0'}  % initial solution point of toolbox
          varargout = [ varargout { efunc.x0(func.x_idx) } ]; %#ok<AGROW>
        case 't0'    % initial tangent of toolbox
          varargout = [ varargout { efunc.tx(func.x_idx) } ]; %#ok<AGROW>
        case {'v0' 'vecs'}  % initial solution point of toolbox - not in use according to Harry
          varargout = [ varargout { efunc.V0(func.x_idx,:) } ]; %#ok<AGROW>
        otherwise
          error('%s: unknown function data field ''%s''', ...
            mfilename, varargin{oarg});
      end
    end
  else
    for oarg = 1:nargin-2
      switch lower(varargin{oarg})
        case {'x0' 'u0'} % full solution vector constructed so far
          varargout = [ varargout { efunc.x0 } ]; %#ok<AGROW>
        case {'xidx' 'uidx'} % indices of full solution vector constructed so far
          varargout = [ varargout { 1:numel(efunc.x0) } ]; %#ok<AGROW>
        case 't0' % full initial tangent vector constructed so far
          varargout = [ varargout { efunc.tx } ]; %#ok<AGROW>
        case 'pidx' % indices of active continuation parameters
          assert(isfield(efunc, 'p_idx'), ...
            '%s: cannot extract %s, equations not closed', ...
            varargin{oarg}, mfilename);
          varargout = [ varargout { efunc.p_idx(:) } ]; %#ok<AGROW>
        case 'fidx' % list of positions of functions added so far
          list = {};
          for i=1:numel(efunc.funcs)
            func = efunc.funcs(i);
            if ~isempty(func.f_idx)
              list = [ list ; { func.identifyer func.f_idx(:) } ]; %#ok<AGROW>
            end
          end
          varargout = [ varargout { list } ]; %#ok<AGROW>
        case {'x0idx' 'u0idx'} % indices of added part of solution vector
          list = {};
          xidx = [];
          for i=1:numel(efunc.funcs)
            func  = efunc.funcs(i);
            x0idx = setdiff(func.x_idx, xidx);
            if ~isempty(x0idx)
              list = [ list ; { func.identifyer x0idx(:) } ]; %#ok<AGROW>
            end
            xidx  = union(xidx, x0idx);
          end
          varargout = [ varargout { list } ]; %#ok<AGROW>
        otherwise
          error('%s: unknown function data field ''%s''', ...
            mfilename, varargin{oarg});
      end
    end
  end
end
