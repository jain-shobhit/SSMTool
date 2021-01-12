function prob = coco_add_glue(prob, fid, x1idx, x2idx, varargin)
%COCO_ADD_GLUE   Add gluing conditions.
%
% PROB = COCO_ADD_GLUE(PROB, FID, X1_IDX, X2_IDX, VARARGIN)
%
% VARARGIN = { [GAP] [PAR_NAMES [PAR_TYPE]] }
%
% PROB      - Continuation problem structure
% FID       - Function identifier
% X1_IDX    - Index array
% X2_IDX    - Index array
% GAP       - Scalar
% PAR_NAMES - String or cell array of strings that label corresponding
%             continuation parameters.
% PAR_TYPE  - One of 'zero', 'active', 'inactive', 'internal', 'regular',
%             or 'singular. Defaults to 'inactive'.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_glue.m 2967 2017-01-10 18:51:26Z hdankowicz $

assert(numel(x1idx)==numel(x2idx), ...
  '%s: number of elements in x1 and x2 must be equal', mfilename);

par_names = {};
par_type  = 'zero';
gap       = zeros(numel(x1idx),1);
argidx    = 1;
if nargin>=argidx+4 && isnumeric(varargin{argidx})
  gap    = varargin{argidx};
  argidx = argidx + 1;
end
if nargin>=argidx+4
  par_names = varargin{argidx};
  par_type  = 'inactive';
  argidx    = argidx + 1;
end
if nargin>=argidx+4
  par_type  = varargin{argidx};
end

data.x1idx = 1:numel(x1idx);
data.x2idx = numel(x1idx) + (1:numel(x2idx));
data.gap   = gap;
data.J     = [speye(numel(x1idx)) , -speye(numel(x2idx))];
xidx       = [ x1idx(:) ; x2idx(:) ];

switch par_type
  case 'zero'
    prob = coco_add_func(prob, fid, @func_GLUE, @func_DGLUEDX, ...
      @func_DGLUEDXDX, data, 'zero', 'xidx', xidx);
  otherwise
    prob = coco_add_func(prob, fid, @func_GLUE, @func_DGLUEDX, ...
      @func_DGLUEDXDX, data, par_type, par_names, 'xidx', xidx);
end

end

function [data, f] = func_GLUE(prob, data, u) %#ok<INUSL>
f = data.gap + data.J*u;
end

function [data, J] = func_DGLUEDX(prob, data, u)  %#ok<INUSD,INUSL>
J = data.J;
end

function [data, dJ] = func_DGLUEDXDX(prob, data, u)  %#ok<INUSD,INUSL>

[m,n] = size(data.J);
dJ = zeros(m,n,n);

end
