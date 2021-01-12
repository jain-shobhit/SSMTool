function hspo_arg_check(tbid, fhan, dfdxhan, dfdphan, modes, events, ...
  resets, t0, x0, p0, pnames)
%HSPO_ARG_CHECK   Basic argument checking for 'hspo' toolbox.
%
% Validate user-supplied inputs and terminate execution with suitable error
% messages if the inputs fail to be of the correct type.
%
% Identical to hspo_v1.
%
% HSPO_ARG_CHECK(TBID, FHAN, DFDXHAN, DFDPHAN, MODES, EVENTS, ...
%   RESETS, T0, X0, P0, PNAMES)
%
% TBID    - Toolbox instance identifier.
% FHAN    - Cell array of function handles to vector field, event function,
%           and reset function.
% DFDXHAN - Optional cell array of function handles to Jacobians w.r.t.
%           problem variables.
% DFDPHAN - Optional cell array of function handles to Jacobians w.r.t.
%           problem parameters.
% MODES   - Sequence of mode identifiers (cell).
% EVENTS  - Sequence of event identifiers (cell).
% RESETS  - Sequence of reset identifiers (cell).
% T0      - Collection of arrays of temporal mesh points (cell).
% X0      - Collection of array of state vectors at mesh points.
% P0      - Initial solution guess for parameter values.
% PNAMES  - String label or cell array of string labels for
%           continuation parameters tracking problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: hspo_arg_check.m 2839 2015-03-05 17:09:01Z fschild $

assert(numel(fhan)==3, ...
  '%s: incomplete input for ''fhan''', tbid);
assert(all(cellfun('isclass', fhan, 'function_handle')), ...
  '%s: input for ''fhan'' not function handles', tbid);
assert(numel(dfdxhan)==3, ...
  '%s: incomplete input for ''dfdxhan''', tbid);
assert(numel(dfdphan)==3, ...
  '%s: incomplete input for ''dfdphan''', tbid);
nos = [numel(modes) numel(events) numel(resets) numel(t0) numel(x0)];
assert(~any(diff(nos)), ...
  '%s incompatible segment specification', tbid);
assert(numel(p0)==numel(pnames) || isempty(pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  tbid);

end
