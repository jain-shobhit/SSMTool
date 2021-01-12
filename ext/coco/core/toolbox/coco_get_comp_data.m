function varargout = coco_get_comp_data(prob, fid, varargin)

% Copyright (C) Harry Dankowicz
% $Id: coco_get_comp_data.m 2973 2017-01-10 19:48:17Z hdankowicz $

complementary = prob.complementary;
varargout = {};
eidx      = find(strcmpi(fid, efunc.identifyers), 1);
oidx      = find(strcmpi(fid, ['complementary' complementary.identifyers]), 1)-1;
if isempty([eidx oidx])
  error('%s: could not find function with identifyer ''%s''', mfilename, fid);
elseif ~isempty(eidx)
  func = coco_get_func_data(prob, fid, 'data');
  for oarg = 1:nargin-2
    switch lower(varargin{oarg})
      case 'data'  % function data structure
        varargout = [ varargout { func.data } ]; %#ok<AGROW>
      case 'uidx'  % K_u function dependency index set
        varargout = [ varargout { func.u_idx(:) } ]; %#ok<AGROW>
      case 'lidx'  % K_l function dependency index set
        varargout = [ varargout { func.l_idx(:) } ]; %#ok<AGROW>
      case 'vidx'  % K_v function dependency index set
        varargout = [ varargout { func.u_idx(:) } ]; %#ok<AGROW>
      case 'type'  % function type
        varargout = [ varargout { func.type } ]; %#ok<AGROW>
      case 'fidx'  % indices to individual functions, empty if non-embedded
        varargout = [ varargout { func.f_idx(:) } ]; %#ok<AGROW>
      case 'v0'    % initial values of complementary continuation variables
        varargout = [ varargout { complementary.v0(func.v_idx) } ]; %#ok<AGROW>
      case 'tv0'    % initial rates of change of complementary continuation variables
        varargout = [ varargout { complementary.tv(func.v_idx) } ]; %#ok<AGROW>
      otherwise
        error('%s: unknown function data field ''%s''', ...
          mfilename, varargin{oarg});
    end
  end
else
  if oidx
    func = complementary.funcs(oidx);
    for oarg = 1:nargin-2
      switch lower(varargin{oarg})
        case 'data'  % function data structure
          if isa(func.data, 'coco_func_data')
            varargout = [ varargout { func.data.protect() } ]; %#ok<AGROW>
          else
            varargout = [ varargout { func.data } ]; %#ok<AGROW>
          end
        case 'uidx'  % K_u function dependency index set
          varargout = [ varargout { func.u_idx(:) } ]; %#ok<AGROW>
        case 'lidx'  % K_l function dependency index set
          varargout = [ varargout { func.l_idx(:) } ]; %#ok<AGROW>
        case 'vidx'  % K_v function dependency index set
          varargout = [ varargout { func.u_idx(:) } ]; %#ok<AGROW>
        case 'type'  % function type
          varargout = [ varargout { func.type } ]; %#ok<AGROW>
        case 'fidx'  % indices to individual functions, empty if non-embedded
          varargout = [ varargout { func.f_idx(:) } ]; %#ok<AGROW>
        case 'v0'    % initial values of complementary continuation variables
          varargout = [ varargout { complementary.v0(func.v_idx) } ]; %#ok<AGROW>
        case 'tv0'    % initial rates of change of complementary continuation variables
          varargout = [ varargout { complementary.tv(func.v_idx) } ]; %#ok<AGROW>
        otherwise
          error('%s: unknown function data field ''%s''', ...
            mfilename, varargin{oarg});
      end
    end
  else
    for oarg = 1:nargin-2
      switch lower(varargin{oarg})
        case 'v0'   % full solution vector constructed so far
          varargout = [ varargout { complementary.v0 } ]; %#ok<AGROW>
        case 'vidx' % indices of full solution vector constructed so far
          varargout = [ varargout { 1:complementary.vdim } ]; %#ok<AGROW>
        case 'tv0' % full initial tangent vector constructed so far
          varargout = [ varargout { complementary.tv } ]; %#ok<AGROW>
        case 'fidx' % list of positions of functions added so far
          list = {};
          for i=1:numel(complementary.funcs)
            func = complementary.funcs(i);
            if ~isempty(func.f_idx)
              list = [ list ; { func.identifyer func.f_idx(:) } ]; %#ok<AGROW>
            end
          end
          varargout = [ varargout { list } ]; %#ok<AGROW>
        case 'v0idx' % indices of added part of solution vector
          list = {};
          vidx = [];
          for i=1:numel(complementary.funcs)
            func  = complementary.funcs(i);
            v0idx = setdiff(func.v_idx, vidx);
            if ~isempty(v0idx)
              list = [ list ; { func.identifyer v0idx(:) } ]; %#ok<AGROW>
            end
            vidx  = union(vidx, v0idx);
          end
          varargout = [ varargout { list } ]; %#ok<AGROW>
        otherwise
          error('%s: unknown function data field ''%s''', ...
            mfilename, varargin{oarg});
      end
    end
  end
end
