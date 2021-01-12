function prob = coco_add_comp_pars(prob, fid, pidx, par_names, varargin)
%COCO_ADD_COMP_PARS   Add external parameters to continuation problem.
%
% COCO_ADD_COMP_PARS(PROB, FID, PIDX, PAR_NAMES VARARGIN)
%
% VARARGIN = { [PAR_TYPE='inactive'] }
% Assign parameter names to continuation multipliers.

% Copyright (C) Harry Dankowicz
% $Id: coco_add_pars.m 2971 2017-01-10 19:41:31Z hdankowicz $

if isempty(fid)
  fid = 'comp_pars';
end

if nargin<5
  par_type = 'inactive';
else
  par_type = varargin{1};
end

assert(~strcmpi(par_type, 'zero'), ...
  '%s: parameter type must be one of active, inactive, internal, regular or singular', ...
  mfilename);

prob = coco_add_comp(prob, fid, @func_PARS, @func_DPARSDX, [], ...
  par_type, par_names, 'lidx', pidx);

end

function [data, y] = func_PARS(~, data, ~, l, ~)
y = l;
end

function [data, J] = func_DPARSDX(~, data, ~, l, ~)
J = {[], speye(size(l,1)), []};
end
