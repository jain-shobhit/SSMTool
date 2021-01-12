function alg_arg_check(data, x0, p0)
%ALG_ARG_CHECK   Basic argument checking for 'alg' toolbox.
%
% Validate user-supplied inputs and terminate execution with suitable error
% messages if the inputs fail to be of the correct type.
%
% ALG_ARG_CHECK(DATA, X0, P0)
%
% DATA - Toolbox data structure.
% X0   - Initial solution guess for problem variables.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_arg_check.m 2839 2015-03-05 17:09:01Z fschild $

assert(isa(data.fhan, 'function_handle'), ...
  'alg: input for ''f'' is not a function handle');
assert(isnumeric(x0), 'alg: input for ''x0'' is not numeric');
assert(isnumeric(p0), 'alg: input for ''p0'' is not numeric');
assert(numel(p0)==numel(data.pnames) || isempty(data.pnames), ...
  'alg: incompatible number of elements for ''p0'' and ''pnames''');

end
