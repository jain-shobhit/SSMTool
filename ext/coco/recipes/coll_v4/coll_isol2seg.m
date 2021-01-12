function prob = coll_isol2seg(prob, oid, varargin)
%COLL_ISOL2SEG   Append 'coll' instance constructed from initial data.
%
% Parse input sequence to construct toolbox data and initial solution guess
% and use this to construct an instance of 'coll'.
%
% Differs from coll_v1 by providing support for modified call to
% coll_init_data.
%
% PROB     = COLL_ISOL2SEG(PROB, OID, VARARGIN)
% VARARGIN = { @F [(@DFDX | '[]') [(@DFDP | '[]')]] T0 X0 [PNAMES] P0 }
%
% PROB   - Continuation problem structure.
% OID    - Object instance identifier (string).
% @F     - Function handle to vector field.
% @DFDX  - Optional function handle to Jacobian w.r.t. problem variables.
% @DFDP  - Optional function handle to Jacobian w.r.t. problem parameters.
% T0     - Array of temporal mesh points.
% X0     - Array of state vectors at mesh points (rectangular matrix with
%          size(x0,1) == numel(t0).
% PNAMES - Optional string label or cell array of string labels for
%          continuation parameters tracking problem parameters.
% P0     - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_isol2seg.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'coll'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
data.fhan = str.get;
data.dfdxhan  = [];
data.dfdphan  = [];
if is_empty_or_func(str.peek)
  data.dfdxhan = str.get;
  if is_empty_or_func(str.peek)
    data.dfdphan = str.get;
  end
end
t0 = str.get;
x0 = str.get;
data.pnames = {};
if iscellstr(str.peek('cell'))
  data.pnames = str.get('cell');
end
p0 = str.get;

coll_arg_check(tbid, data, t0, x0, p0);     % Validate input
data = coll_get_settings(prob, tbid, data); % Get toolbox settings
data = coll_init_data(data, t0, x0, p0);    % Build toolbox data
sol  = coll_init_sol(data, t0, x0, p0);     % Build initial solution guess
prob = coll_construct_seg(prob, tbid, data, sol); % Append continuation problem

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end
