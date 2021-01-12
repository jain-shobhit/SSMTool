function prob = dft_sol2orb(prob, oid, varargin)
%DFT_SOL2SEG   Append 'dft' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% PROB     = DFT_SOL2SEG(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB ['end_dft']}
%
% PROB - Continuation problem structure.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: dft_sol2orb.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'dft'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run  = str.get;
if ischar(str.peek)
  soid = str.get;
  lab  = str.get;
else
  soid = oid;
  lab  = str.get;
end
if strcmp(str.peek, 'end-dft')
  str.skip;
end

[sol data] = dft_read_solution(soid, run, lab);  % Extract solution and toolbox data from disk
data       = dft_get_settings(prob, tbid, data); % Get toolbox settings
sol        = dft_init_sol(data, sol.t, sol.x, sol.p);  % Build initial solution guess
data       = dft_init_data(data, sol);                 % Build toolbox data
prob       = dft_construct_orb(prob, tbid, data, sol); % Append continuation problem

end
