function prob = coco_add_pars(prob, fid, varargin)
%COCO_ADD_PARS   Add external parameters to continuation problem.
%
% COCO_ADD_PARS(PROB, FID, VARARGIN)
%
% VARARGIN = { PIDX PAR_NAMES [PAR_TYPE='inactive'] }
% Assign parameter names to variables.
%
% VARARGIN = { PAR_NAMES PVALS [TVALS] [PAR_TYPE='inactive'] }
% Add variables as parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_pars.m 2971 2017-01-10 19:41:31Z hdankowicz $

if isempty(fid)
  fid = 'pars';
end

if isnumeric(varargin{1})
  if nargin<5
    varargin{3} = 'inactive';
  end
  prob = map_vars(prob, fid, varargin{:});
else
  if nargin<5
    varargin{3} = [];
    varargin{4} = 'inactive';
  elseif nargin<6
    if ischar(varargin{3})
      varargin{4} = varargin{3};
      varargin{3} = [];
    else
      varargin{4} = 'inactive';
    end
  end
  prob = add_pars(prob, fid, varargin{:});
end

end

function prob = map_vars(prob, fid, pidx, par_names, par_type)

assert(~strcmpi(par_type, 'zero'), ...
  '%s: parameter type must be one of active, inactive, internal, regular or singular', ...
  mfilename);

if ischar(par_names)
  par_names = { par_names };
end

prob = coco_add_func(prob, fid, @func_PARS, @func_DPARSDX, ...
  @func_DPARSDXDX, [], par_type, par_names, 'xidx', pidx);

end

function prob = add_pars(prob, fid, par_names, pvals, tvals, par_type)

assert(~any(strcmpi(par_type, {'zero' 'regular' 'singular'})), ...
  '%s: parameter type must be one of active, inactive or internal', ...
  mfilename);

if ischar(par_names)
  par_names = { par_names };
end

prob = coco_add_func(prob, fid, @func_PARS, @func_DPARSDX, ...
  @func_DPARSDXDX, [], par_type, par_names, 'x0', pvals, 't0', tvals);

end

function [data, g] = func_PARS(prob, data, xp) %#ok<INUSL>
g = xp;
end

function [data, J] = func_DPARSDX(prob, data, xp)  %#ok<INUSL>
J = speye(size(xp,1));
end

function [data, dJ] = func_DPARSDXDX(prob, data, xp)  %#ok<INUSL>

nx = size(xp,1);
dJ = zeros(nx,nx,nx);

end
