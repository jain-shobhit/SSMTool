function prob = alg_isol2eqn(prob, oid, varargin)
%ALG_ISOL2EQN   Append 'alg' instance constructed from initial data.
%
% Parse input sequence to construct toolbox data and initial solution guess
% and use this to construct an instance of 'alg'.
%
% Differs from alg_v5 by providing support for optional settings governing
% the execution of 'alg' constructors and function objects.
%
% PROB     = ALG_ISOL2EQN(PROB, OID, VARARGIN)
% VARARGIN = { @F [(@DFDX |'[]') [(@DFDP |'[]')]] X0 [PNAMES] P0 }
%
% PROB   - Continuation problem structure.
% OID    - Object instance identifier (string).
% @F     - Function handle to zero function.
% @DFDX  - Optional function handle to Jacobian w.r.t. problem variables.
% @DFDP  - Optional function handle to Jacobian w.r.t. problem parameters.
% X0     - Initial solution guess for problem variables.
% PNAMES - Optional string label or cell array of string labels for
%          continuation parameters tracking problem parameters.
% P0     - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_isol2eqn.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'alg');  % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
data.fhan = str.get;
data.dfdxhan = [];
data.dfdphan = [];
if is_empty_or_func(str.peek)
  data.dfdxhan = str.get;
  if is_empty_or_func(str.peek)
    data.dfdphan = str.get;
  end
end
x0 = str.get;
data.pnames = {};
if iscellstr(str.peek('cell'))
  data.pnames = str.get('cell');
end
p0 = str.get;

alg_arg_check(tbid, data, x0, p0);          % Validate input
data  = alg_get_settings(prob, tbid, data); % Get toolbox settings
data  = alg_init_data(data, x0, p0);        % Build toolbox data
sol.u = [x0(:); p0(:)];                     % Build initial solution guess
prob  = alg_construct_eqn(prob, tbid, data, sol); % Append continuation problem

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end
