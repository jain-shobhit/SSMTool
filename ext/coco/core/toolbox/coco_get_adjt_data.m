function varargout = coco_get_adjt_data(prob, fid, varargin)

% Copyright (C) Harry Dankowicz
% $Id: coco_get_adjt_data.m 2972 2017-01-10 19:44:34Z hdankowicz $

adjoint   = prob.adjoint;
varargout = {};
idx       = find(strcmpi(fid, ['adjoint' adjoint.identifyers]), 1)-1;
if isempty(idx)
  error('%s: could not find adjoint with identifyer ''%s''', mfilename, fid);
else
  if idx
    func = adjoint.funcs(idx);
    for oarg = 1:nargin-2
      switch lower(varargin{oarg})
        case 'data'  % adjoint function data structure
          if isa(func.data, 'coco_func_data')
            varargout = [ varargout { func.data.protect() } ]; %#ok<AGROW>
          else
            varargout = [ varargout { func.data } ]; %#ok<AGROW>
          end
        case 'afidx' % rows of adjoint function matrix associated with adjoint function
          varargout = [ varargout { func.af_idx(:) } ]; %#ok<AGROW>
        case 'axidx' % columns of adjoint function matrix associated with adjoint function
          varargout = [ varargout { func.ax_idx(:) } ]; %#ok<AGROW>
        case 'l0'    % initial values for continuation multipliers associated with adjoint function
          varargout = [ varargout { adjoint.l0(func.af_idx) } ]; %#ok<AGROW>
        case 'tl0'    % initial values for rates of change of continuation multipliers associated with adjoint function
          varargout = [ varargout { adjoint.tl0(func.af_idx) } ]; %#ok<AGROW>
        otherwise
          error('%s: unknown function data field ''%s''', ...
            mfilename, varargin{oarg});
      end
    end
  else
    for oarg = 1:nargin-2
      switch lower(varargin{oarg})
        case 'l0' % initial values for all continuation multipliers
          varargout = [ varargout { adjoint.l0 } ]; %#ok<AGROW>
        case 'lidx' % indices of all continuation multipliers
          varargout = [ varargout { 1:numel(adjoint.l0) } ]; %#ok<AGROW>
        case 'lidx_phi' % indices of continuation multipliers associated with zero functions
          varargout = [ varargout { setdiff(1:numel(adjoint.l0), ...
            [adjoint.l_idx adjoint.s_idx]) } ]; %#ok<AGROW>
        case 'lidx_psi' % indices of continuation multipliers associated with inactive or active monitor functions
          varargout = [ varargout { adjoint.l_idx } ]; %#ok<AGROW>
        case 'lidx_g' % indices of continuation multipliers associated with inactive or active monitor functions
          varargout = [ varargout { adjoint.s_idx } ]; %#ok<AGROW>
        case 'tl0' % initial rates of change for all continuation multipliers
          varargout = [ varargout { adjoint.tl0 } ]; %#ok<AGROW>
        case 'afidx' % indices of rows of adjoint matrix associated with individual adjoint functions
          list = {};
          for i=1:numel(adjoint.funcs)
            func = adjoint.funcs(i);
            if ~isempty(func.af_idx)
              list = [ list ; { func.identifyer func.af_idx(:) } ]; %#ok<AGROW>
            end
          end
          varargout = [ varargout { list } ]; %#ok<AGROW>
        case 'axidx' % indices of added columns of adjoint matrix associated with individual adjoint functions
          list = {};
          xidx = [];
          for i=1:numel(adjoint.funcs)
            func  = adjoint.funcs(i);
            axidx = setdiff(func.ax_idx, xidx);
            if ~isempty(axidx)
              list = [ list ; { func.identifyer axidx(:) } ]; %#ok<AGROW>
            end
            xidx  = union(xidx, axidx);
          end
          varargout = [ varargout { list } ]; %#ok<AGROW>
        otherwise
          error('%s: unknown function data field ''%s''', ...
            mfilename, varargin{oarg});
      end
    end
  end
end
