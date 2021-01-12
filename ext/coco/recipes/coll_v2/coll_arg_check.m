function coll_arg_check(tbid, data, t0, x0, p0)
%COLL_ARG_CHECK   Basic argument checking for 'coll' toolbox.
%
% Validate user-supplied inputs and terminate execution with suitable error
% messages if the inputs fail to be of the correct type.
%
% Identical to coll_v1.
%
% COLL_ARG_CHECK(TBID, DATA, T0, X0, P0)
%
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% T0   - Array of temporal mesh points.
% X0   - Array of state vectors at mesh points.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_arg_check.m 2839 2015-03-05 17:09:01Z fschild $

assert(isa(data.fhan, 'function_handle'), ...
  '%s: input for ''f'' is not a function handle', tbid);
assert(isnumeric(t0), '%s: input for ''t0'' is not numeric', tbid);
assert(isnumeric(x0), '%s: input for ''x0'' is not numeric', tbid);
assert(isnumeric(p0), '%s: input for ''p0'' is not numeric', tbid);
assert(ndims(t0)==2 && min(size(t0))==1, ...
  '%s: input for ''t0'' is not a vector', tbid);
assert(ndims(x0)==2, ...
  '%s: input for ''x0'' is not an array of vectors', tbid);
assert(size(x0, 1)==numel(t0), ...
  '%s: dimensions of ''t0'' and ''x0'' do not match', tbid);
assert(numel(p0)==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  tbid);

end
