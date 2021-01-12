function prob = alg_sol2eqn(prob, oid, varargin)
%ALG_SOL2EQN   Append 'alg' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% Differs from alg_v4 by providing support for modifying an existing
% continuation problem structure, and relying on unique object instance
% identifiers to distinguish between toolbox instances.
%
% PROB     = ALG_SOL2EQN(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB     - Continuation problem structure.
% OID      - Target object instance identifier (string).
% RUN      - Run identifier (string).
% SOID     - Source object instance identifier (string).
% LAB      - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_sol2eqn.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'alg');  % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run  = str.get;
if ischar(str.peek)
  soid = str.get;
else
  soid = oid;
end
lab = str.get;

[sol data] = alg_read_solution(soid, run, lab); % Extract solution and toolbox data from disk
prob       = alg_construct_eqn(prob, tbid, data, sol); % Append continuation problem

end
