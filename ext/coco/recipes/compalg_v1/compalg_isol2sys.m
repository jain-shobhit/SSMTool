function prob = compalg_isol2sys(prob, oid, varargin)
%COMPALG_ISOL2SYS   Append and glue several 'alg' instances constructed from initial data.
%
% Use multiple calls to alg_isol2eqn and close the system to ensure that
% problem parameters are identical across all instances of 'alg'.
%
% PROB     = COMPALG_ISOL2SYS(PROB, OID, VARARGIN)
% VARARGIN = { FUNCS {X0...} [PNAMES] P0 }
% FUNCS    = {@F...} [ {(@DFDX | '[]')...} [ {(@DFDP | '[]')...} ] ]
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
% $Id: compalg_isol2sys.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'compalg'); % Create toolbox instance identifier
str  = coco_stream(varargin{:});    % Convert varargin to stream of tokens for argument parsing
fhans = str.get('cell');            % Parse cell array or single token into cell array
data.neqs = numel(fhans);           % Number of 'alg' objects
dfdxhans = cell(1, data.neqs);
dfdphans = cell(1, data.neqs);
if is_empty_or_func(str.peek('cell'))
  dfdxhans = str.get('cell');       % Parse cell array or single token into cell array
  if is_empty_or_func(str.peek('cell'))
    dfdphans = str.get('cell');     % Parse cell array or single token into cell array
  end
end
x0 = str.get('cell');               % Parse cell array or single token into cell array
data.pnames = {};
if iscellstr(str.peek('cell'))
  data.pnames = str.get('cell');
end
p0 = str.get;

compalg_arg_check(tbid, data, dfdxhans, dfdphans, x0, p0); % Validate input
for i=1:data.neqs
  toid = coco_get_id(tbid, sprintf('eqn%d', i)); % Create unique 'alg' object instance identifier
  prob = alg_isol2eqn(prob, toid, fhans{i}, dfdxhans{i}, ...
    dfdphans{i}, x0{i}, p0);                     % Call 'alg' constructor without pnames
end
prob = compalg_close_sys(prob, tbid, data);      % Append gluing conditions

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = all(cellfun('isempty', x) | cellfun('isclass', x, 'function_handle'));
end
