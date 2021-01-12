function compalg_arg_check(tbid, data, dfdxhans, dfdphans, x0, p0)
%COMPALG_ARG_CHECK   Basic argument checking for 'compalg' toolbox.
%
% Validate user-supplied inputs and terminate execution with suitable error
% messages if the inputs fail to be of the correct type.
%
% COMPALG_ARG_CHECK(TBID, DATA, DFDXHANS, DFDPHANS, X0, P0)
%
% TBID     - Toolbox instance identifier.
% DATA     - Toolbox data structure.
% DFDXHANS - (cell array of) Function handles to Jacobian(s) of zero
%            function(s) w.r.t. problem variables.
% DFDPHANS - (cell array of) Function handles to Jacobian(s) of zero
%            function(s) w.r.t. problem parameters.
% X0       - (cell array of) Initial solution guess(es) for problem
%            variables.
% P0       - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: compalg_arg_check.m 2839 2015-03-05 17:09:01Z fschild $

assert(data.neqs~=0, '%s: insufficient number of equations', tbid);
assert(data.neqs==numel(dfdxhans), ...
  '%s: incompatible number of inputs in ''dfdxhans''', tbid);
assert(data.neqs==numel(dfdphans), ...
  '%s: incompatible number of inputs in ''dfdphans''', tbid);
assert(data.neqs==numel(x0), ...
  '%s: incompatible number of inputs in ''x0''', tbid);
assert(numel(p0)==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  tbid);

end
