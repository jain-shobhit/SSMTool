function tor_arg_check(tbid, data)
%COLL_ARG_CHECK   Basic argument checking for 'tor' toolbox.
%
% Validate user-supplied inputs and terminate execution with suitable error
% messages if the inputs fail to be of the correct type.
%
% TOR_ARG_CHECK(TBID, DATA)
%
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.

assert(isa(data.fhan, 'function_handle'), ...
  '%s: input for ''f'' is not a function handle', tbid);
assert(is_empty_or_func(data.dfdxhan), ...
    '%s: input for ''fx'' is neither empty nor a function handle', tbid);
assert(is_empty_or_func(data.dfdphan), ...
    '%s: input for ''fp'' is neither empty nor a function handle', tbid);
assert(is_empty_or_func(data.dfdthan), ...
    '%s: input for ''ft'' is neither empty nor a function handle', tbid);
assert(isnumeric(data.t0), '%s: input for ''T0'' is not numeric', tbid);
assert(isnumeric(data.x0), '%s: input for ''X0''  is not numeric', tbid);
assert(isnumeric(data.p0), '%s: input for ''p0'' is not numeric', tbid);
assert(numel(data.p0)==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  tbid);
end


function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end
