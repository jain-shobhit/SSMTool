function [prob, data] = coco_add_functionals(prob, fid, A, b, x_idx, varargin)
%COCO_ADD_FUNCTIONALS   Add linear functionals to zero problem.
%
% [PROB [DATA]] = COCO_ADD_FUNCTIONALS(PROB, FID, A, B, X_IDX, VARARGIN)
%
% VARARGIN = { [PAR_NAMES, TYPE] }
%
% 1) A (m,n)-matrix, X_IDX (n)-vector   -> A*x=B
% 2) A (m,n)-matrix, X_IDX (m,n)-matrix -> sum(A.*x,2)=B

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_functionals.m 2969 2017-01-10 19:13:51Z hdankowicz $

%% parse input arguments
%  varargin = { type pnames }

if nargin>5
  par_names = varargin{1};
  par_type  = varargin{2};
  if ischar(par_names)
    par_names = { par_names };
  elseif isnumeric(par_names)
    par_names = coco_get_def_par_names('PAR', par_names);
  end
else
  par_type  = 'zero';
  par_names = [];
end

[m, n]      = size(A);
data.A      = A;
data.b      = b;
data.xshape = [n m];
data.cols   = prod(data.xshape);
N = numel(x_idx);
if m*n==N % sum(A.*x,2)=B
  data.xidx = reshape(1:N, m, n);
elseif n==N % A*x=B
  data.xidx = repmat(1:N, m, 1);
  x_idx     = repmat(x_idx(:)', m, 1);
else
  error('%s: inconsistent dimensions', mfilename);
end
data.fidx   = repmat( (1:m)', 1, n);

if nargout>1
  data = coco_func_data(data);
end

switch par_type
  case 'zero'
    prob = coco_add_func(prob, fid, @func_LFUNC, @func_DLFUNCDX, ...
      @func_DLFUNCDXDX, data, 'zero', 'xidx', x_idx');
  otherwise
    prob = coco_add_func(prob, fid, @func_LFUNC, @func_DLFUNCDX, ...
      @func_DLFUNCDXDX, data, par_type, par_names, 'xidx', x_idx');
end

end

function [data, g] = func_LFUNC(prob, data, xp) %#ok<INUSL>
%FUNC_LFUNC  Evaluate linear functional(s).

x = reshape(xp, data.xshape)';
g = sum(data.A.*x, 2)-data.b;

end

function [data, J] = func_DLFUNCDX(prob, data, xp)  %#ok<INUSD,INUSL>
%FUNC_DLFUNCDX  Compute linearisation of linear functional(s).
J = sparse(data.fidx, data.xidx, data.A, data.xshape(2), data.cols);
end

function [data, dJ] = func_DLFUNCDXDX(prob, data, xp) %#ok<INUSD,INUSL>
%FUNC_DLFUNCDXDX  Compute second derivatives of linear functional(s).
dJ = zeros(data.xshape(2), data.cols, data.cols);
end
