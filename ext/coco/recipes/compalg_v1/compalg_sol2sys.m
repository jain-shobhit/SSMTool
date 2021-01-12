function prob = compalg_sol2sys(prob, oid, varargin)
%COMPALG_SOL2SYS   Append and glue several 'alg' instances constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% PROB     = COMPALG_SOL2SYS(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB - Continuation problem structure.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: compalg_sol2sys.m 2839 2015-03-05 17:09:01Z fschild $

ttbid = coco_get_id(oid, 'compalg'); % Create toolbox instance identifier
str = coco_stream(varargin{:});      % Convert varargin to stream of tokens for argument parsing
run = str.get;
if ischar(str.peek)
  stbid = coco_get_id(str.get, 'compalg');
else
  stbid = ttbid;
end
lab = str.get;

data = coco_read_solution(stbid, run, lab);       % Extract 'compalg' toolbox data from disk
for i=1:data.neqs
  soid = coco_get_id(stbid, sprintf('eqn%d', i)); % Create 'alg' source object instance identifier
  toid = coco_get_id(ttbid, sprintf('eqn%d', i)); % Create 'alg' target object instance identifier
  prob = alg_sol2eqn(prob, toid, run, soid, lab); % Call 'alg' constructor
end
prob = compalg_close_sys(prob, ttbid, data);      % Append gluing conditions

end
