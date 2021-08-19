function prob = ode_tor2tor(prob, oid, varargin)
%ODE_tor2tor   Append 'tor' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% PROB     = ODE_TOR2TOR(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB - Continuation problem structure.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

tbid = coco_get_id(oid, 'tor'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run  = str.get;
if ischar(str.peek)
  soid = str.get;
else
  soid = tbid;
end
lab = str.get;

prob = ode_bvp2bvp(prob, tbid, run, soid, lab); % Append continuation problem
end
