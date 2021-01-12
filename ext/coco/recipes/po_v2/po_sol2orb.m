function prob = po_sol2orb(prob, oid, varargin)
%PO_SOL2ORB   Append 'po' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% Differs from po_v1 by supporting optional toolbox settings.
%
% PROB     = PO_SOL2ORB(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB - Continuation problem structure.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_sol2orb.m 2839 2015-03-05 17:09:01Z fschild $

ttbid = coco_get_id(oid, 'po');   % Create toolbox instance identifier
str   = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run   = str.get;
if ischar(str.peek)
  stbid = coco_get_id(str.get, 'po');
else
  stbid = ttbid;
end
lab = str.get;

toid = coco_get_id(ttbid, 'seg'); % Target 'coll' object instance identifier
soid = coco_get_id(stbid, 'seg'); % Source 'coll' object instance identifier
prob = coll_sol2seg(prob, toid, run, soid, lab); % Construct 'coll' instance
data = coco_read_solution(stbid, run, lab);
data = po_get_settings(prob, tbid, data); % Get toolbox settings
data = po_init_data(prob, ttbid, data);   % Build toolbox data
prob = po_close_orb(prob, ttbid, data);   % Append boundary conditions

end
