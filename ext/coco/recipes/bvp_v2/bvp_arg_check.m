function bvp_arg_check(tbid, data)
%BVP_ARG_CHECK   Basic argument checking for 'bvp' toolbox.
%
% Validate user-supplied inputs and terminate execution with suitable error
% messages if the inputs fail to be of the correct type.
%
% Identical to bvp_v1.
%
% BVP_ARG_CHECK(TBID, DATA)
%
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: bvp_arg_check.m 2839 2015-03-05 17:09:01Z fschild $

assert(isa(data.fhan, 'function_handle'), ...
  '%s: input for ''f'' is not a function handle', tbid);

end
