function prob = compalg_isol2sys(prob, oid, varargin)
%COMPALG_ISOL2SYS   Append and glue several 'alg' instances constructed from initial data.
%
% Use multiple calls to alg_isol2eqn and close the system to ensure that
% problem parameters are identical across all instances of 'alg'.
%
% Differs from compalg_v1 by use of recursive argument parsing.
%
% PROB     = COMPALG_ISOL2SYS(PROB, OID, VARARGIN)
% VARARGIN = { ALG... [ PNAMES | 'end-alg' ] }
%
% PROB   - Continuation problem structure.
% OID    - Object instance identifier (string).
% ALG    - Argument to alg_isol2eqn without PNAMES
% PNAMES - Optional string label or cell array of string labels for
%          continuation parameters tracking problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: compalg_isol2sys.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'compalg'); % Create toolbox instance identifier
str  = coco_stream(varargin{:});    % Convert varargin to stream of tokens for argument parsing
data.neqs = 0;
while isa(str.peek, 'function_handle')
  data.neqs = data.neqs+1;
  toid      = coco_get_id(tbid, sprintf('eqn%d', data.neqs)); % Create unique object instance identifier
  prob      = alg_isol2eqn(prob, toid, str); % Use 'alg' constructor to parse one instance of ALG
end
data.pnames = {};
if strcmpi(str.peek, 'end-alg') % Check for stop token
  str.skip;
elseif iscellstr(str.peek('cell'))
  data.pnames = str.get('cell');
end

compalg_arg_check(prob, tbid, data);        % Validate input
prob = compalg_close_sys(prob, tbid, data); % Append gluing conditions

end
