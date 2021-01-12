function alg_arg_check(tbid, data, x0, p0)
%ALG_ARG_CHECK   Basic argument checking for 'alg' toolbox.
%
% Validate user-supplied inputs and terminate execution with suitable error
% messages if the inputs fail to be of the correct type.
%
% Identical to alg_v5.
%
% ALG_ARG_CHECK(TBID, DATA, X0, P0)
%
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% X0   - Initial solution guess for problem variables.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_arg_check.m 2839 2015-03-05 17:09:01Z fschild $

assert(isa(data.fhan, 'function_handle'), ...
  '%s: input for ''f'' is not a function handle', tbid);
assert(isnumeric(x0), '%s: input for ''x0'' is not numeric', tbid);
assert(isnumeric(p0), '%s: input for ''p0'' is not numeric', tbid);
assert(numel(p0)==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  tbid);

end
